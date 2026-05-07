# -*- coding: utf-8 -*-
"""
Script 01 — Generate route point shapefiles.

Reads the freight flow OD matrix and the road network, then for every
long-haul OD pair (straight-line distance > MIN_HAUL_KM) generates a
sequence of points spaced POINT_SPACING_KM kilometres apart along the
routed road path. Results are written as one shapefile per OD pair to
ROUTE_POINTS_DIR.

Points carry the attributes of the nearest road edge (speed class,
traffic counts, edge ID) which are used in the next pipeline step to
apply rest/break rules and compute charging demand.

If the script is interrupted it can be safely restarted: existing output
shapefiles are detected and skipped automatically.
"""

import math
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Geod, Transformer
from shapely.geometry import LineString, Point
from shapely import ops

# Allow running this script from any working directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    FLOW_CSV, NUTS3_CSV, NODES_CSV, EDGES_CSV, CENTROIDS_SHP,
    ROUTE_POINTS_DIR, POINT_SPACING_KM, MIN_HAUL_KM,
    MIN_FLOW_TONS_2030, N_WORKERS,
)

# World Mercator projection used for the straight-line distance filter.
# This projection exaggerates distances at high latitudes, but using it
# here ensures the same set of routes is selected as in the original model.
_to_mercator = Transformer.from_crs('EPSG:4326', 'EPSG:3395', always_xy=True)

# Geodesic object for accurate point spacing along WGS84 ellipsoid
_geod = Geod(ellps='WGS84')

pd.options.mode.chained_assignment = None

# ---------------------------------------------------------------------------
# Worker-process globals
# ---------------------------------------------------------------------------
# Each worker process receives these indexed DataFrames once at pool startup
# via _init_worker. Storing them as module globals avoids re-pickling large
# DataFrames on every individual task call.
nodes_idx     = None
edges_idx     = None
centroids_idx = None


def _init_worker(n_idx, e_idx, c_idx):
    global nodes_idx, edges_idx, centroids_idx
    nodes_idx     = n_idx
    edges_idx     = e_idx
    centroids_idx = c_idx


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mercator_distance_km(lon1, lat1, lon2, lat2):
    x1, y1 = _to_mercator.transform(lon1, lat1)
    x2, y2 = _to_mercator.transform(lon2, lat2)
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) / 1000


def _as_scalar(series_or_scalar):
    if isinstance(series_or_scalar, pd.Series):
        return series_or_scalar.iloc[0]
    return series_or_scalar


def _build_edge_line(edge_id):
    """
    Look up a network edge by ID and return (attribute_dict, LineString).
    Returns None if the edge or either endpoint node is missing.
    """
    try:
        er = edges_idx.loc[edge_id]
        if isinstance(er, pd.DataFrame):
            er = er.iloc[0]
        sn = nodes_idx.loc[int(er['Network_Node_A_ID'])]
        en = nodes_idx.loc[int(er['Network_Node_B_ID'])]
        if isinstance(sn, pd.DataFrame): sn = sn.iloc[0]
        if isinstance(en, pd.DataFrame): en = en.iloc[0]
    except KeyError:
        return None

    geom = LineString([
        (sn['Network_Node_X'], sn['Network_Node_Y']),
        (en['Network_Node_X'], en['Network_Node_Y']),
    ])
    attr = {'Network_Edge_ID': edge_id}
    attr.update(er.to_dict())
    return attr, geom


def _points_along_geodesic(geom, spacing_m):
    """
    Walk a LineString in EPSG:4326 and return Point objects spaced every
    spacing_m geodesic metres along the WGS84 ellipsoid.

    Geodesic spacing is used because projected (Mercator) interpolation
    inflates distances at European latitudes by up to 56% at 50°N, which
    would place points only ~6.4 km apart instead of 10 km.
    """
    from shapely.geometry import MultiLineString
    if isinstance(geom, MultiLineString):
        open_parts = [p for p in geom.geoms
                      if list(p.coords)[0] != list(p.coords)[-1]]
        parts = open_parts if open_parts else list(geom.geoms)
        return _points_along_geodesic(max(parts, key=lambda p: p.length), spacing_m)

    coords = list(geom.coords)
    if not coords:
        return []
    if len(coords) == 1:
        return [Point(coords[0])]

    result    = [Point(coords[0])]
    remainder = 0.0

    for i in range(len(coords) - 1):
        lon1, lat1 = coords[i]
        lon2, lat2 = coords[i + 1]
        az12, _, seg_len = _geod.inv(lon1, lat1, lon2, lat2)
        offset = spacing_m - remainder
        while offset <= seg_len:
            lon, lat, _ = _geod.fwd(lon1, lat1, az12, offset)
            result.append(Point(lon, lat))
            offset += spacing_m
        remainder = seg_len - (offset - spacing_m)

    return result


def _generate_points(merged_gdf):
    spacing_m = POINT_SPACING_KM * 1000
    points = []
    for geom in merged_gdf.geometry:
        points.extend(_points_along_geodesic(geom, spacing_m))
    if not points:
        return None
    return gpd.GeoDataFrame(geometry=points, crs='EPSG:4326')


def _resolve_multiline(multi, origin_centroid, dest_centroid):
    """
    Reconstruct a single origin-to-destination path from a MultiLineString.

    linemerge produces a MultiLineString when edges form a star topology —
    multiple segments meeting at a shared junction node. This function
    identifies the junction, separates origin-side from destination-side
    segments, and assembles the correct through-path.
    """
    from collections import Counter

    origin_pt = Point(origin_centroid['Geometric_center_X'],
                      origin_centroid['Geometric_center_Y'])
    dest_pt   = Point(dest_centroid['Geometric_center_X'],
                      dest_centroid['Geometric_center_Y'])

    open_parts = [p for p in multi.geoms
                  if list(p.coords)[0] != list(p.coords)[-1]]
    if not open_parts:
        open_parts = list(multi.geoms)

    ep_count: Counter = Counter()
    for p in open_parts:
        ep_count[p.coords[0]]  += 1
        ep_count[p.coords[-1]] += 1
    junction, junction_freq = ep_count.most_common(1)[0]

    if junction_freq < 2:
        # No shared junction — pick the part whose endpoints span origin→dest
        best, best_score = None, float('inf')
        for part in open_parts:
            first, last = Point(part.coords[0]), Point(part.coords[-1])
            score = min(origin_pt.distance(first) + dest_pt.distance(last),
                        origin_pt.distance(last)  + dest_pt.distance(first))
            if score < best_score:
                best_score, best = score, part
        return best

    origin_side, dest_side = [], []
    for part in open_parts:
        far_coords = (list(part.coords) if part.coords[0] == junction
                      else list(reversed(list(part.coords))))
        far_end = Point(far_coords[-1])
        if origin_pt.distance(far_end) < dest_pt.distance(far_end):
            origin_side.append((origin_pt.distance(far_end), far_coords))
        else:
            dest_side.append((dest_pt.distance(far_end), far_coords))

    if origin_side and dest_side:
        origin_side.sort(key=lambda x: x[0])
        dest_side.sort(key=lambda x: x[0])
        origin_coords = list(reversed(origin_side[0][1]))
        dest_coords   = dest_side[0][1]
        return LineString(origin_coords + dest_coords[1:])

    all_coords = [x[1] for x in (origin_side or dest_side)]
    all_coords.sort(key=lambda c: origin_pt.distance(Point(c[-1])))
    return LineString(list(reversed(all_coords[0])))


def _orient_origin_to_dest(points_gdf, origin_id):
    """
    Reverse the point sequence if it runs destination→origin.

    Network edges have no guaranteed traversal direction, so the merged
    geometry may be oriented backwards. Comparing the first and last point
    against the origin centroid detects this and corrects it.
    """
    try:
        cr = centroids_idx.loc[int(origin_id)]
        if isinstance(cr, pd.DataFrame): cr = cr.iloc[0]
    except KeyError:
        return points_gdf

    origin_pt = Point(cr['Geometric_center_X'], cr['Geometric_center_Y'])
    if (origin_pt.distance(points_gdf.iloc[-1].geometry)
            < origin_pt.distance(points_gdf.iloc[0].geometry)):
        points_gdf = points_gdf.iloc[::-1].reset_index(drop=True)
    return points_gdf


# ---------------------------------------------------------------------------
# Per-route worker function
# ---------------------------------------------------------------------------

def process_route(row):
    """
    Build and save the point shapefile for one OD pair.

    Steps:
      1. Skip if the output file already exists.
      2. Skip if the straight-line distance is below MIN_HAUL_KM.
      3. Build LineString geometries for each edge in the routed path.
      4. Merge all edge segments into one continuous LineString.
      5. Generate evenly spaced points along the merged line.
      6. Orient the sequence so it runs origin → destination.
      7. Assign the nearest edge attributes to each point via spatial join.
      8. Write the result as a shapefile.
    """
    origin = int(row['ID_origin_region'])
    dest   = int(row['ID_destination_region'])
    pair   = '{}_{}'.format(origin, dest)
    out_path = os.path.join(ROUTE_POINTS_DIR, 'route_points_{}.shp'.format(pair))

    if os.path.exists(out_path):
        return

    try:
        sc = centroids_idx.loc[origin]
        ec = centroids_idx.loc[dest]
        if isinstance(sc, pd.DataFrame): sc = sc.iloc[0]
        if isinstance(ec, pd.DataFrame): ec = ec.iloc[0]
    except KeyError:
        print('missing centroid: {}'.format(pair))
        return

    dist_km = _mercator_distance_km(
        sc['Geometric_center_X'], sc['Geometric_center_Y'],
        ec['Geometric_center_X'], ec['Geometric_center_Y'],
    )
    if dist_km <= MIN_HAUL_KM:
        return

    edges_str = row['Edge_path_E_road']
    if not isinstance(edges_str, str):
        try:
            edges_str = edges_str.values[0]
        except Exception:
            edges_str = str(edges_str)

    records = []
    try:
        for edge_name in edges_str[1:-1].split(','):
            edge_name = edge_name.strip()
            if not edge_name:
                continue
            result = _build_edge_line(int(edge_name))
            if result is not None:
                attr, geom = result
                attr['geometry'] = geom
                records.append(attr)

        if not records:
            return

        lines_gdf = gpd.GeoDataFrame(records, crs='EPSG:4326')
        lines_gdf.rename(columns={
            'Network_Edge_ID': 'Edge_ID',
            'Network_Node_A_ID': 'A_ID',
            'Network_Node_B_ID': 'B_ID',
        }, inplace=True)

        union_geom  = lines_gdf.geometry.union_all()
        merged_geom = (union_geom if union_geom.geom_type == 'LineString'
                       else ops.linemerge(union_geom))

        if merged_geom.geom_type == 'MultiLineString':
            merged_geom = _resolve_multiline(merged_geom, sc, ec)

        merged_gdf  = gpd.GeoDataFrame(geometry=[merged_geom], crs='EPSG:4326')
        points_gdf  = _generate_points(merged_gdf)
        if points_gdf is None or points_gdf.empty:
            return

        points_gdf  = _orient_origin_to_dest(points_gdf, origin)

        # Assign road attributes to each point using the nearest edge.
        # 'nearest' is more robust than 'intersects' because floating-point
        # rounding after CRS conversion can shift points fractionally off
        # the exact edge geometry.
        result_gdf = gpd.sjoin_nearest(points_gdf, lines_gdf, how='left')
        result_gdf.to_file(out_path)
        print(pair)

    except Exception as e:
        print('error {}: {}'.format(pair, e))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    # The __name__ guard is required on Windows: multiprocessing uses the
    # "spawn" start method, which re-imports this module in each worker.
    # Without the guard the worker would re-execute the pool creation,
    # triggering an infinite process spawn loop.
    t0 = time.time()
    os.makedirs(ROUTE_POINTS_DIR, exist_ok=True)

    flow = pd.read_csv(FLOW_CSV)
    flow = flow[flow['Traffic_flow_tons_2030'] > MIN_FLOW_TONS_2030]

    _nodes        = pd.read_csv(NODES_CSV)
    nodes_idx     = _nodes.set_index('Network_Node_ID')

    _edges = pd.read_csv(EDGES_CSV)
    _edges.rename(columns={
        'Traffic_flow_trucks_2019': 'Trucks_19',
        'Traffic_flow_trucks_2030': 'Trucks_30',
    }, inplace=True)
    _edges.drop(
        ['Unnamed: 0', 'Unnamed: 0.1', 'Manually_Added', 'Distance'],
        axis=1, inplace=True, errors='ignore',
    )
    edges_idx = _edges.set_index('Network_Edge_ID')

    _centroids    = pd.read_csv(NUTS3_CSV)
    centroids_idx = _centroids.set_index('ETISPlus_Zone_ID')

    rows      = [row for _, row in flow.iterrows()]
    n_workers = N_WORKERS or mp.cpu_count()
    print('Processing {:,} routes across {} workers...'.format(len(rows), n_workers))

    with ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=_init_worker,
        initargs=(nodes_idx, edges_idx, centroids_idx),
    ) as executor:
        list(executor.map(process_route, rows, chunksize=50))

    print('Done in {:.2f} minutes'.format((time.time() - t0) / 60))
