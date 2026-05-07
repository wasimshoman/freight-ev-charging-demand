# -*- coding: utf-8 -*-
"""
Script 02 — Apply EU rest/break rules and compute EV charging demand.

Reads the per-route point shapefiles produced by script 01, then for each
route:
  - Computes point-to-point distances and travel times from road speed data.
  - Flags ferry segments and no-charge edges where charging is unavailable.
  - Walks the point sequence and applies EU Regulation 561/2006 mandatory
    break and rest rules to identify where stops occur.
  - Records the coordinates, energy consumed, and distance driven at each stop.

All stop and rest events are collected into three output shapefiles:
  stop_rest_events.shp  — every stop/rest event across all routes
  breaks.shp            — mandatory 45-min break events only
  rests.shp             — mandatory 9-h rest events only

Each event row carries the OD flow data (annual trucks and tonnes) so that
downstream scripts can weight demand by freight volume.
"""

import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    FLOW_CSV, NUTS3_CSV, CENTROIDS_SHP,
    ROUTE_POINTS_DIR, BASE_OUTPUT_DIR, STOP_REST_FILE,
    MIN_FLOW_TRUCKS_2019, MIN_DISTANCE_KM, MAX_SEGMENT_KM,
    FERRY_SPEED_CODE, BREAK_THRESHOLD_MIN, BREAK_DURATION_MIN,
    REST_DURATION_MIN, ONE_DRIVER_MAX_MIN,
    SINGLE_DRIVER_BREAKS, TWO_DRIVER_BREAKS,
    ENERGY_RATE_LOW_KWH_KM, ENERGY_RATE_HIGH_KWH_KM,
    ENERGY_SPEED_THRESHOLD_KMH, NO_CHARGE_EDGE_IDS,
    N_WORKERS, ANNUAL_WORKING_DAYS, TRUCK_CAPACITY_T,
)

# Column names for the stop/rest event DataFrame
_EVENT_COLUMNS = [
    'ID_origin_region', 'ID_destination_region',
    'YTruck(30)', 'YFlow(30)', 'YTruck(19)', 'YFlow(19)',
    'Name_origin_region', 'Name_destination_region',
    'X', 'Y', 'Rest', 'Break', 'ChaDisKM', 'ChaEnekWh',
]

# ---------------------------------------------------------------------------
# Worker-process globals
# ---------------------------------------------------------------------------
centroids_idx = None


def _init_worker(c_idx):
    global centroids_idx
    centroids_idx = c_idx


# ---------------------------------------------------------------------------
# Type-coercion helpers
# ---------------------------------------------------------------------------

def _to_int(value):
    if isinstance(value, (np.int64, int)):
        return value
    try:
        return value.values[0]
    except Exception:
        return int(value)


def _to_float(value):
    if isinstance(value, (np.float64, float)):
        return value
    try:
        return value.values[0]
    except Exception:
        return float(value)


def _to_str(value):
    if isinstance(value, str):
        return value
    try:
        return value.values[0]
    except Exception:
        return str(value)


def _flow_values(file_row):
    """Extract the eight flow-file fields stored in every stop/rest event."""
    return [
        _to_int(file_row['ID_origin_region']),
        _to_int(file_row['ID_destination_region']),
        _to_float(file_row['Traffic_flow_trucks_2030']),
        _to_float(file_row['Traffic_flow_tons_2030']),
        _to_float(file_row['Traffic_flow_trucks_2019']),
        _to_float(file_row['Traffic_flow_tons_2019']),
        _to_str(file_row['Name_origin_region']),
        _to_str(file_row['Name_destination_region']),
    ]


# ---------------------------------------------------------------------------
# Distance and travel time
# ---------------------------------------------------------------------------

def _add_distances_and_time(gdf, origin_id):
    """
    Compute point-to-point distances (km) and travel times (min) for a route.

    Points are projected to EPSG:3395 for distance calculation to match the
    Mercator-based distances used in the original freight model. The first
    point's distance is measured from the origin zone centroid rather than
    from a preceding point, since the centroid represents the trip start.

    Segments longer than MAX_SEGMENT_KM are geometry gaps from disconnected
    network edges and are zeroed so they do not inflate travel time. Ferry
    segments retain their travel time (the ferry crossing counts toward trip
    duration) but have their road distance set to zero.
    """
    try:
        cr = centroids_idx.loc[int(origin_id)]
        if isinstance(cr, pd.DataFrame):
            cr = cr.iloc[0]
        origin_pt = Point(cr['Geometric_center_X'], cr['Geometric_center_Y'])
    except KeyError:
        origin_pt = gdf.geometry.iloc[0]

    gdf_proj    = gdf.to_crs('EPSG:3395')
    origin_proj = (
        gpd.GeoDataFrame({'geometry': [origin_pt]}, crs='EPSG:4326')
        .to_crs('EPSG:3395')
        .geometry.iloc[0]
    )

    prev_geom       = gdf_proj.geometry.shift()
    prev_geom.iloc[0] = origin_proj

    gdf = gdf.copy()
    gdf['TraDisKm'] = gdf_proj.geometry.distance(prev_geom) / 1000
    gdf.loc[gdf['TraDisKm'] > MAX_SEGMENT_KM, 'TraDisKm'] = 0.0
    gdf['TraTimMin'] = gdf['TraDisKm'] * 60 / gdf['SpeedF'].astype(float)
    gdf.loc[gdf['SpeedF'] == FERRY_SPEED_CODE, 'TraDisKm'] = 0.0

    return gdf


def _mark_no_charge_edges(gdf):
    """Set NoCh=1 for points on edges where charging is not available."""
    gdf = gdf.copy()
    if 'Edge_ID' in gdf.columns:
        gdf['NoCh'] = gdf['Edge_ID'].astype(str).isin(NO_CHARGE_EDGE_IDS).astype(int)
    else:
        gdf['NoCh'] = 0
    return gdf


# ---------------------------------------------------------------------------
# EU rest/break rule engine
# ---------------------------------------------------------------------------

def _apply_rest_break_rules(gdf, reg_travel_time, file_row):
    """
    Walk the point sequence and apply EU Regulation 561/2006 stop rules.

    The estimated trip duration (reg_travel_time) determines the regime:

      Single-driver  (BREAK_THRESHOLD_MIN ≤ trip < ONE_DRIVER_MAX_MIN):
        - SINGLE_DRIVER_BREAKS mandatory breaks of BREAK_DURATION_MIN each
        - followed by a REST_DURATION_MIN daily rest

      Two-driver  (trip ≥ ONE_DRIVER_MAX_MIN):
        - TWO_DRIVER_BREAKS mandatory breaks of BREAK_DURATION_MIN each
        - followed by a REST_DURATION_MIN daily rest

    Ferry segments (SpeedF == FERRY_SPEED_CODE) and no-charge edges
    (NoCh == 1) are passive — they do not accumulate driving time toward
    the next break threshold.

    Returns (events, updated_gdf) where events is a list of stop/rest rows
    for this route and updated_gdf has Stop, Rest, AccTrTi, CharDis,
    CharEner columns added.

    The inner loop extracts all required columns to numpy arrays before
    iterating, and writes results back in bulk after the loop. This avoids
    per-row pandas index lookups which are ~100× slower than numpy indexing
    for sequences of this length.
    """
    is_single_driver = BREAK_THRESHOLD_MIN <= reg_travel_time < ONE_DRIVER_MAX_MIN
    is_two_drivers   = reg_travel_time >= ONE_DRIVER_MAX_MIN

    if not (is_single_driver or is_two_drivers):
        return [], gdf

    n = len(gdf)

    speeds     = gdf['SpeedF'].astype(int).to_numpy()
    trav_times = gdf['TraTimMin'].to_numpy(dtype=float)
    trav_dists = gdf['TraDisKm'].to_numpy(dtype=float)
    noch       = (gdf['NoCh'].to_numpy(dtype=int)
                  if 'NoCh' in gdf.columns else np.zeros(n, dtype=int))
    xs         = gdf.geometry.x.to_numpy()
    ys         = gdf.geometry.y.to_numpy()

    char_dis_arr  = trav_dists.copy()
    char_ener_arr = np.where(
        speeds < ENERGY_SPEED_THRESHOLD_KMH,
        trav_dists * ENERGY_RATE_LOW_KWH_KM,
        trav_dists * ENERGY_RATE_HIGH_KWH_KM,
    )

    stops_arr  = np.zeros(n, dtype=int)
    rests_arr  = np.zeros(n, dtype=int)
    acc_ti_arr = np.zeros(n, dtype=float)

    accumulated_time = 0.0
    trip_time        = 0.0
    pauses           = 0
    used_energy      = 0.0
    char_dist        = 0.0
    events           = []
    fv               = _flow_values(file_row)

    for i in range(n):
        used_energy      += char_ener_arr[i]
        char_dist        += trav_dists[i]
        accumulated_time += trav_times[i]
        trip_time        += trav_times[i]

        stop = 0
        rest = 0

        is_passive = (speeds[i] == FERRY_SPEED_CODE) or (noch[i] == 1)

        if not is_passive:
            if is_single_driver:
                if trip_time > BREAK_THRESHOLD_MIN and pauses < SINGLE_DRIVER_BREAKS:
                    trip_time = 0; pauses += 1; stop = 1
                    accumulated_time += BREAK_DURATION_MIN
                    events.append(fv + [xs[i], ys[i], 0, 1, char_dist, used_energy])
                    used_energy = char_dist = 0.0
                elif trip_time > BREAK_THRESHOLD_MIN and pauses >= SINGLE_DRIVER_BREAKS:
                    trip_time = 0; pauses = 0; stop = 1; rest = 1
                    accumulated_time += REST_DURATION_MIN
                    events.append(fv + [xs[i], ys[i], 1, 0, char_dist, used_energy])
                    used_energy = char_dist = 0.0

            elif is_two_drivers:
                if trip_time > BREAK_THRESHOLD_MIN and pauses < TWO_DRIVER_BREAKS:
                    trip_time = 0; pauses += 1; stop = 1
                    accumulated_time += BREAK_DURATION_MIN
                    events.append(fv + [xs[i], ys[i], 0, 1, char_dist, used_energy])
                    used_energy = char_dist = 0.0
                elif trip_time > BREAK_THRESHOLD_MIN and pauses >= TWO_DRIVER_BREAKS:
                    trip_time = 0; pauses = 0; stop = 1; rest = 1
                    accumulated_time += REST_DURATION_MIN
                    events.append(fv + [xs[i], ys[i], 1, 0, char_dist, used_energy])
                    used_energy = char_dist = 0.0

        stops_arr[i]  = stop
        rests_arr[i]  = rest
        acc_ti_arr[i] = accumulated_time

    gdf = gdf.copy()
    gdf['Stop']     = stops_arr
    gdf['Rest']     = rests_arr
    gdf['AccTrTi']  = acc_ti_arr
    gdf['CharDis']  = char_dis_arr
    gdf['CharEner'] = char_ener_arr

    return events, gdf


# ---------------------------------------------------------------------------
# Per-route worker function
# ---------------------------------------------------------------------------

def _process_route(file_row):
    """
    Process one OD flow record: enrich the route shapefile with distance,
    time, and stop/rest columns, then return stop/rest event rows.
    """
    P1 = _to_int(file_row['ID_origin_region'])
    P2 = _to_int(file_row['ID_destination_region'])
    route_path = os.path.join(ROUTE_POINTS_DIR, 'route_points_{}_{}.shp'.format(P1, P2))

    if not os.path.exists(route_path):
        return []

    reg_travel_time = int(file_row['Total_distance'] / 70) * 60

    try:
        gdf = gpd.read_file(route_path)

        needs_distance = (
            'TraTimMin' not in gdf.columns
            or 'TraDisKm' not in gdf.columns
            or gdf['TraTimMin'].sum() == 0
            or ('TraDisKm' in gdf.columns and (gdf['TraDisKm'] > MAX_SEGMENT_KM).any())
        )
        if needs_distance:
            gdf = _add_distances_and_time(gdf, P1)

        if 'NoCh' not in gdf.columns:
            gdf = _mark_no_charge_edges(gdf)

        if 'CharDis'  not in gdf.columns: gdf['CharDis']  = 0.0
        if 'CharEner' not in gdf.columns: gdf['CharEner'] = 0.0

        events, gdf = _apply_rest_break_rules(gdf, reg_travel_time, file_row)
        gdf.to_file(route_path)
        return events

    except Exception as e:
        print('error {}_{}: {}'.format(P1, P2, e))
        return []


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _add_dtn_columns(gdf):
    """
    Add daily truck-need (DTN) columns at every 10% EV penetration step
    for both 2019 and 2030 annual flows.

    DTN = annual_flow_tonnes * penetration_share / (working_days * payload_per_truck)
    """
    for flow_year, flow_col in [('30', 'YFlow(30)'), ('19', 'YFlow(19)')]:
        for pct in range(10, 110, 10):
            gdf['DTN{}_{}'.format(flow_year, pct)] = round(
                gdf[flow_col] * (pct / 100) / (ANNUAL_WORKING_DAYS * TRUCK_CAPACITY_T), 2
            )
    return gdf


def _build_output_gdf(events_df, ref_crs):
    events_df = events_df.astype({
        'YTruck(30)': int, 'YFlow(30)': int,
        'YTruck(19)': int, 'YFlow(19)': int,
    })
    gdf = gpd.GeoDataFrame(
        events_df,
        geometry=[Point(x, y) for x, y in zip(events_df['X'], events_df['Y'])],
        crs=ref_crs,
    )
    return _add_dtn_columns(gdf)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    t0 = time.time()
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    flow = pd.read_csv(FLOW_CSV)
    flow = flow[flow['Traffic_flow_trucks_2019'] > MIN_FLOW_TRUCKS_2019]
    flow = flow[flow['Total_distance'] > MIN_DISTANCE_KM]

    _centroids    = pd.read_csv(NUTS3_CSV)
    centroids_idx = _centroids.set_index('ETISPlus_Zone_ID')

    rows      = [row for _, row in flow.iterrows()]
    n_workers = N_WORKERS or mp.cpu_count()
    print('Processing {:,} routes across {} workers...'.format(len(rows), n_workers))

    with ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=_init_worker,
        initargs=(centroids_idx,),
    ) as executor:
        results = list(executor.map(_process_route, rows, chunksize=50))

    all_events = [ev for route_events in results for ev in route_events]

    if not all_events:
        print('No stop/rest events found — verify that route shapefiles exist in ROUTE_POINTS_DIR.')
    else:
        out_df  = pd.DataFrame(all_events, columns=_EVENT_COLUMNS)
        ref_crs = gpd.read_file(CENTROIDS_SHP).crs

        out_dir = os.path.dirname(STOP_REST_FILE)
        os.makedirs(out_dir, exist_ok=True)

        _build_output_gdf(out_df, ref_crs).to_file(STOP_REST_FILE)
        _build_output_gdf(out_df[out_df['Break'] == 1].copy(), ref_crs).to_file(
            os.path.join(out_dir, 'breaks.shp'))
        _build_output_gdf(out_df[out_df['Rest'] == 1].copy(), ref_crs).to_file(
            os.path.join(out_dir, 'rests.shp'))

        print('Saved stop_rest_events.shp, breaks.shp, rests.shp to:')
        print(' ', out_dir)

    print('Done in {:.2f} minutes'.format((time.time() - t0) / 60))
