# -*- coding: utf-8 -*-
"""
Script 05 — Aggregate truck charging demand by origin zone and country.

Re-derives daily truck demand directly from the annual OD flow file
(rather than from the pre-aggregated stop event data) and aggregates it
first to origin zone level, then to country level via spatial join.

Outputs:
  origin_truck_demand.shp   — one point per origin zone
    TMainDTN  daily trucks needing a charge (main EV scenario)
    TDTN100   daily trucks needing a charge (100% electrification)

  country_truck_demand.gpkg — one polygon per country
    TMainDTN  summed across all origin zones within the country
    TDTN100   summed across all origin zones within the country

The distinction from script 04 is in how demand is computed:
  Script 04 sums the DTN columns already embedded in each stop event row.
  Script 05 re-computes DTN from raw annual tonnage, which reflects the
  full flow volume for each OD pair regardless of where stops occurred.
"""

import os
import sys
import time

import geopandas as gpd
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    STOP_REST_FILE, FLOW_CSV, CENTROIDS_SHP, EU_COUNTRIES_SHP,
    BASE_OUTPUT_DIR, MIN_FLOW_TRUCKS_2019, MIN_DISTANCE_KM,
    ANNUAL_WORKING_DAYS, TRUCK_CAPACITY_T, EV_SHARE_MAIN,
)

OUTPUT_POINTS = os.path.join(BASE_OUTPUT_DIR, 'origin_truck_demand.shp')
OUTPUT_COUNTRY = os.path.join(BASE_OUTPUT_DIR, 'country_truck_demand.gpkg')


if __name__ == '__main__':
    t0 = time.time()

    # Step 1 — unique OD pairs present in the stop/rest event data
    print('Loading unique OD pairs from stop/rest events...')
    needed = ['ID_origin_', 'ID_destina']
    try:
        rb = gpd.read_file(STOP_REST_FILE, columns=needed, ignore_geometry=True)
    except TypeError:
        rb = gpd.read_file(STOP_REST_FILE)
        rb = pd.DataFrame(rb[needed])

    rb['ID_origin_'] = rb['ID_origin_'].astype(int)
    rb['ID_destina'] = rb['ID_destina'].astype(int)
    od_pairs = rb.drop_duplicates(subset=['ID_origin_', 'ID_destina'])
    print('  {:,} unique OD pairs from {:,} origin zones'.format(
        len(od_pairs), od_pairs['ID_origin_'].nunique()))

    # Step 2 — annual tonnage per OD pair from the flow file
    print('Loading flow file...')
    flow = pd.read_csv(FLOW_CSV)
    flow = flow[flow['Traffic_flow_trucks_2019'] > MIN_FLOW_TRUCKS_2019]
    flow = flow[flow['Total_distance'] > MIN_DISTANCE_KM]
    flow = flow[['ID_origin_region', 'ID_destination_region',
                 'Traffic_flow_tons_2030']].copy()
    flow = flow.rename(columns={
        'ID_origin_region':       'ID_origin_',
        'ID_destination_region':  'ID_destina',
        'Traffic_flow_tons_2030': 'YearlyTons',
    })
    print('  {:,} flow rows after filtering'.format(len(flow)))

    # Step 3 — compute daily truck need per OD pair
    # DTN = annual_tonnes / working_days / payload_per_truck
    merged = od_pairs.merge(flow, on=['ID_origin_', 'ID_destina'], how='inner')
    merged['MainDTN'] = merged['YearlyTons'] / ANNUAL_WORKING_DAYS * EV_SHARE_MAIN / TRUCK_CAPACITY_T
    merged['DTN100']  = merged['YearlyTons'] / ANNUAL_WORKING_DAYS / TRUCK_CAPACITY_T
    print('  {:,} OD pairs matched'.format(len(merged)))

    # Step 4 — group by origin zone
    grouped = (
        merged.groupby('ID_origin_')
        .agg(TMainDTN=('MainDTN', 'sum'), TDTN100=('DTN100', 'sum'))
        .reset_index()
    )

    # Step 5 — attach centroid geometry
    print('Attaching centroid geometry...')
    centroids = gpd.read_file(CENTROIDS_SHP)[['ETISPlus_Z', 'geometry']]
    centroids = centroids.rename(columns={'ETISPlus_Z': 'ID_origin_'})

    origin_gdf = gpd.GeoDataFrame(
        grouped.merge(centroids, on='ID_origin_', how='inner'),
        geometry='geometry',
        crs=centroids.crs,
    )
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)
    origin_gdf.to_file(OUTPUT_POINTS)
    print('  {:,} origin points saved: {}'.format(
        len(origin_gdf), os.path.basename(OUTPUT_POINTS)))

    # Step 6 — aggregate to country level via spatial join
    print('Aggregating to country level...')
    countries = gpd.read_file(EU_COUNTRIES_SHP).to_crs(origin_gdf.crs)

    joined = gpd.sjoin(
        origin_gdf,
        countries[['name', 'geometry']],
        how='inner',
        predicate='intersects',
    )

    country_agg = (
        joined.groupby('name')
        .agg(TMainDTN=('TMainDTN', 'sum'), TDTN100=('TDTN100', 'sum'))
        .reset_index()
        .rename(columns={'name': 'Country'})
    )

    country_geom = countries[['name', 'geometry']].rename(columns={'name': 'Country'})
    gpd.GeoDataFrame(
        country_agg.merge(country_geom, on='Country', how='left'),
        geometry='geometry',
        crs=countries.crs,
    ).to_file(OUTPUT_COUNTRY, driver='GPKG')
    print('  {:,} countries saved: {}'.format(
        len(country_agg), os.path.basename(OUTPUT_COUNTRY)))

    print('\nDone in {:.2f} minutes'.format((time.time() - t0) / 60))
