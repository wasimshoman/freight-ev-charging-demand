# -*- coding: utf-8 -*-
"""
Script 03 — Aggregate stop/rest events to a 25×25 km grid and size chargers.

Reads the stop/rest events from script 02 and produces spatial summaries
at two levels:

  Grid level (charging_demand_grid.shp):
    Each populated 25×25 km cell reports the estimated number of fast
    chargers (30-min and 60-min CCS) and slow chargers (overnight MCS)
    required to serve the EV truck demand at that location. Charger counts
    are derived from a pre-computed M/G/c queuing model table.

  Country level (country_summary.csv / .gpkg):
    Grid cells are spatially joined to country boundaries and aggregated
    to give a per-country breakdown of charging infrastructure demand and
    daily energy consumption.

The REFERENCE_GRID_SHP setting in config.py controls whether a fixed
reference grid (recommended for reproducibility) or a fresh auto-generated
grid is used for spatial aggregation.
"""

import math
import os
import sys
import time

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    STOP_REST_FILE, BASE_OUTPUT_DIR,
    QUEUE_30MIN_CSV, QUEUE_60MIN_CSV, EU_COUNTRIES_SHP,
    REFERENCE_GRID_SHP, GRID_SIZE_M, GRID_CRS,
    EV_SHARE_MAIN, ANNUAL_WORKING_DAYS, TRUCK_CAPACITY_T,
)


# ---------------------------------------------------------------------------
# Queuing model lookup
# ---------------------------------------------------------------------------

def _lookup_chargers(lambda_val, queue_df):
    """
    Return the minimum number of charger bays for the given hourly arrival
    rate, using the pre-computed M/G/c queuing model table.

    For arrival rates beyond the table range the result is extrapolated in
    chunks of 152 arrivals/hour (the table's maximum modelled rate).
    """
    if lambda_val <= 0:
        return 0
    matches = queue_df[queue_df['lambda'] > lambda_val]
    if not matches.empty:
        return int(matches.iloc[0]['servers'])
    whole     = int(lambda_val // 152)
    remaining = lambda_val % 152
    tail      = queue_df[queue_df['lambda'] > remaining]
    extra     = int(tail.iloc[0]['servers']) if not tail.empty else 0
    return int(80 * whole + extra)


# ---------------------------------------------------------------------------
# Step 1 — enrich stop/rest events with demand columns
# ---------------------------------------------------------------------------

def _enrich_events(gdf):
    """
    Derive charging demand columns from the flow and energy data.

    MainDTN  — daily trucks needing a charge at this stop (main EV scenario)
    ChEMain  — energy charged per stop event (MWh, main scenario)
    ChE30M   — weighted energy (MWh) = ChEMain × MainDTN
    ChERM    — energy at rest stops with slow chargers (MWh)
    ChEBM    — energy at break stops with fast chargers (MWh)
    MDTN_R   — daily trucks using slow chargers (rest stops)
    MDTN_B   — daily trucks using fast chargers (break stops)
    ChER100  — energy at rest stops at 100% electrification (MWh)
    ChEB100  — energy at break stops at 100% electrification (MWh)
    DTN100_R — daily trucks at rest stops, 100% electrification
    DTN100_B — daily trucks at break stops, 100% electrification
    """
    gdf = gdf.copy()
    gdf['MainDTN'] = gdf['YFlow(30)'] / ANNUAL_WORKING_DAYS * EV_SHARE_MAIN / TRUCK_CAPACITY_T
    gdf['ChEMain'] = pd.to_numeric(gdf['ChaEnekWh']) / 1000
    gdf['ChE30M']  = gdf['ChEMain'] * gdf['MainDTN']
    gdf['Rest']    = pd.to_numeric(gdf['Rest'])
    gdf['Break']   = pd.to_numeric(gdf['Break'])
    gdf['ChERM']   = gdf['Rest']  * gdf['ChE30M']
    gdf['MDTN_R']  = gdf['Rest']  * gdf['MainDTN']
    gdf['ChEBM']   = gdf['Break'] * gdf['ChE30M']
    gdf['MDTN_B']  = gdf['Break'] * gdf['MainDTN']
    gdf['ChER100'] = gdf['Rest']  * gdf['ChEMain'] * gdf['DTN30_100']
    gdf['ChEB100'] = gdf['Break'] * gdf['ChEMain'] * gdf['DTN30_100']
    gdf['DTN100_R']= gdf['Rest']  * gdf['YFlow(30)'] / ANNUAL_WORKING_DAYS / TRUCK_CAPACITY_T
    gdf['DTN100_B']= gdf['Break'] * gdf['YFlow(30)'] / ANNUAL_WORKING_DAYS / TRUCK_CAPACITY_T
    return gdf


# ---------------------------------------------------------------------------
# Step 2 — build or load the aggregation grid
# ---------------------------------------------------------------------------

def _make_grid(gdf_proj, cell_size):
    """Generate a regular square grid covering the bounding box of gdf_proj."""
    xmin, ymin, xmax, ymax = gdf_proj.total_bounds
    xmin = math.floor(xmin / cell_size) * cell_size
    ymin = math.floor(ymin / cell_size) * cell_size
    xmax = math.ceil(xmax  / cell_size) * cell_size
    ymax = math.ceil(ymax  / cell_size) * cell_size

    cells = []
    x = xmin
    while x < xmax:
        y = ymin
        while y < ymax:
            cells.append(box(x, y, x + cell_size, y + cell_size))
            y += cell_size
        x += cell_size

    grid = gpd.GeoDataFrame({'geometry': cells}, crs=GRID_CRS)
    grid['cell_id'] = range(len(grid))
    return grid


# ---------------------------------------------------------------------------
# Step 3 — spatial join and aggregate
# ---------------------------------------------------------------------------

def _aggregate_to_grid(rb_gdf, grid_gdf, queue_30, queue_60):
    """
    Spatially join stop/rest events to grid cells and compute per-cell
    aggregates and charger counts.
    """
    rb_proj = rb_gdf.to_crs(GRID_CRS)
    joined  = gpd.sjoin(
        rb_proj,
        grid_gdf[['cell_id', 'geometry']],
        how='left',
        predicate='within',
    )

    agg = (
        joined.dropna(subset=['cell_id'])
        .groupby('cell_id')
        .agg(
            SUM_DTN30_ = ('DTN30_100', 'sum'),
            SUM_DTN301 = ('DTN30_30',  'sum'),
            SUM_ChaEne = ('ChaEnekWh', 'sum'),
            SUM_MainDT = ('MainDTN',   'sum'),
            SUM_ChE30M = ('ChE30M',    'sum'),
            SUM_ChERM  = ('ChERM',     'sum'),
            SUM_MDTN_R = ('MDTN_R',    'sum'),
            SUM_ChEBM  = ('ChEBM',     'sum'),
            SUM_MDTN_B = ('MDTN_B',    'sum'),
            SUM_ChER10 = ('ChER100',   'sum'),
            SUM_ChEB10 = ('ChEB100',   'sum'),
        )
        .round(2)
        .reset_index()
    )

    grid_out = grid_gdf.merge(agg, on='cell_id', how='inner').drop(columns=['cell_id'])

    grid_out['DTN100']  = grid_out['SUM_DTN30_']
    grid_out['DTN30']   = grid_out['SUM_DTN301']
    grid_out['ChE100']  = grid_out['SUM_ChaEne'] / 1000
    grid_out['MainDTN'] = grid_out['SUM_MainDT']
    grid_out['ChE30']   = grid_out['SUM_ChE30M']
    grid_out['ChERM']   = grid_out['SUM_ChERM']
    grid_out['MDTN_R']  = grid_out['SUM_MDTN_R']
    grid_out['ChEBM']   = grid_out['SUM_ChEBM']
    grid_out['MDTN_B']  = grid_out['SUM_MDTN_B']

    # Slow chargers: sized to serve rest-stop trucks at 2 charges per day
    grid_out['NSCh2pD'] = np.ceil(grid_out['MDTN_R'] / 2).astype(int)
    # Fast chargers: arrival rate = trucks × peak-hour share (6%), sized via M/G/c
    grid_out['NFCh30m'] = grid_out['MDTN_B'].apply(
        lambda v: _lookup_chargers(v * 0.06, queue_30))
    grid_out['TotCha']  = grid_out['NSCh2pD'] + grid_out['NFCh30m']
    grid_out['NFCh1H']  = grid_out['MDTN_B'].apply(
        lambda v: _lookup_chargers(v * 0.06, queue_60))
    # 100% electrification scales arrival rate by 1/EV_SHARE_MAIN
    scale_100 = 1.0 / EV_SHARE_MAIN
    grid_out['NFCh100'] = grid_out['MDTN_B'].apply(
        lambda v: _lookup_chargers(v * scale_100 * 0.06, queue_30))
    grid_out['NSCh100'] = np.ceil(grid_out['MDTN_R'] * scale_100 / 2).astype(int)

    return grid_out.set_crs(GRID_CRS)


# ---------------------------------------------------------------------------
# Step 4 — country-level summary
# ---------------------------------------------------------------------------

def _summarise_by_country(grid_gdf, countries_gdf):
    grid_4326 = grid_gdf.to_crs('EPSG:4326')
    countries = countries_gdf.to_crs('EPSG:4326')

    joined = gpd.sjoin(
        grid_4326,
        countries[['name', 'geometry']],
        how='inner',
        predicate='intersects',
    )

    grouped = (
        joined.groupby('name')
        .agg(
            NumberOfStations               = ('geometry', 'count'),
            DailyTrucksAtStops_15pct       = ('MainDTN',  'sum'),
            ChargedEnergy_MWh_100pct       = ('ChE100',   'sum'),
            DailyTrucks_100pct             = ('DTN100',   'sum'),
            ChargedEnergy_MWh_15pct        = ('ChE30',    'sum'),
            SlowChargerEnergy_MWh          = ('ChERM',    'sum'),
            DailyTrucksSlowChargers        = ('MDTN_R',   'sum'),
            FastChargerEnergy_MWh          = ('ChEBM',    'sum'),
            DailyTrucksFastChargers        = ('MDTN_B',   'sum'),
            SlowChargingPoints             = ('NSCh2pD',  'sum'),
            FastChargingPoints_30min       = ('NFCh30m',  'sum'),
            TotalChargingPoints            = ('TotCha',   'sum'),
            FastChargingPoints_60min       = ('NFCh1H',   'sum'),
            FastChargingPoints_30min_100pct= ('NFCh100',  'sum'),
            SlowChargingPoints_100pct      = ('NSCh100',  'sum'),
        )
        .reset_index()
        .rename(columns={'name': 'Country'})
    )

    safe_div = lambda a, b: (a / b.replace(0, np.nan)).round(2)
    grouped['TrucksPerFastCharger'] = safe_div(grouped['DailyTrucksFastChargers'],
                                               grouped['FastChargingPoints_30min'])
    grouped['TrucksPerSlowCharger'] = safe_div(grouped['DailyTrucksSlowChargers'],
                                               grouped['SlowChargingPoints'])
    grouped['SlowToFastRatio']      = safe_div(grouped['SlowChargingPoints'],
                                               grouped['FastChargingPoints_30min'])
    return grouped, joined


def _print_summary(grid_gdf, grouped):
    print('\n=== Main scenario ({:.0f}% electrification) ==='.format(EV_SHARE_MAIN * 100))
    print('  Fast chargers (30 min) : {:>8,.0f}'.format(grid_gdf['NFCh30m'].sum()))
    print('  Slow chargers (2/day)  : {:>8,.0f}'.format(grid_gdf['NSCh2pD'].sum()))
    print('  Total chargers         : {:>8,.0f}'.format(grid_gdf['TotCha'].sum()))
    print('  Fast chargers (1 hr)   : {:>8,.0f}'.format(grid_gdf['NFCh1H'].sum()))
    print('\n=== 100% electrification ===')
    print('  Fast chargers (30 min) : {:>8,.0f}'.format(grid_gdf['NFCh100'].sum()))
    print('  Slow chargers          : {:>8,.0f}'.format(grid_gdf['NSCh100'].sum()))

    fast = grouped['FastChargingPoints_30min'].sum()
    slow = grouped['SlowChargingPoints'].sum()
    n    = grouped['NumberOfStations'].sum()
    print('\n  Avg fast chargers per station : {:.2f}'.format(fast / n))
    print('  Avg slow chargers per station : {:.2f}'.format(slow / n))

    top5 = ['Germany', 'France', 'Poland', 'Spain', 'Italy']
    top5_e = grouped[grouped['Country'].isin(top5)]['ChargedEnergy_MWh_15pct'].sum()
    total_e = grouped['ChargedEnergy_MWh_15pct'].sum()
    print('\n  Top-5 country energy share      : {:.0f}%'.format(top5_e / total_e * 100))
    print('  Total daily energy (GWh, {:.0f}%) : {:.1f}'.format(
        EV_SHARE_MAIN * 100, total_e / 1000))
    print('  Total daily energy (GWh, 100%)  : {:.1f}'.format(
        grouped['ChargedEnergy_MWh_100pct'].sum() / 1000))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    t0 = time.time()
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    print('Loading inputs...')
    rb_gdf    = gpd.read_file(STOP_REST_FILE)
    queue_30  = pd.read_csv(QUEUE_30MIN_CSV)
    queue_60  = pd.read_csv(QUEUE_60MIN_CSV)
    countries = gpd.read_file(EU_COUNTRIES_SHP)
    print('  {:,} stop/rest events loaded'.format(len(rb_gdf)))

    print('\nStep 1 - enriching events with demand columns...')
    rb_gdf = _enrich_events(rb_gdf)
    enriched_path = os.path.join(BASE_OUTPUT_DIR, 'stop_rest_events_enriched.shp')
    rb_gdf.to_file(enriched_path)
    print('  Saved: {}'.format(os.path.basename(enriched_path)))

    if REFERENCE_GRID_SHP and os.path.exists(REFERENCE_GRID_SHP):
        print('\nStep 2 - loading reference grid ({})...'.format(
            os.path.basename(REFERENCE_GRID_SHP)))
        grid_gdf = gpd.read_file(REFERENCE_GRID_SHP).to_crs(GRID_CRS)
        grid_gdf = grid_gdf[['geometry']].copy()
        grid_gdf['cell_id'] = range(len(grid_gdf))
        print('  {:,} cells loaded'.format(len(grid_gdf)))
    else:
        print('\nStep 2 - generating {:,} m grid...'.format(GRID_SIZE_M))
        rb_proj  = rb_gdf.to_crs(GRID_CRS)
        grid_gdf = _make_grid(rb_proj, GRID_SIZE_M)
        print('  {:,} candidate cells'.format(len(grid_gdf)))

    print('\nStep 3 - spatial join and aggregation...')
    agg_gdf  = _aggregate_to_grid(rb_gdf, grid_gdf, queue_30, queue_60)
    grid_path = os.path.join(BASE_OUTPUT_DIR, 'charging_demand_grid.shp')
    agg_gdf.to_crs('EPSG:4326').to_file(grid_path)
    print('  {:,} populated cells saved: {}'.format(
        len(agg_gdf), os.path.basename(grid_path)))

    print('\nStep 4 - country-level summary...')
    grouped, _ = _summarise_by_country(agg_gdf, countries)

    csv_path = os.path.join(BASE_OUTPUT_DIR, 'country_summary.csv')
    grouped.to_csv(csv_path, index=False)
    print('  Saved: {}'.format(os.path.basename(csv_path)))

    countries_4326 = (
        countries.to_crs('EPSG:4326')[['name', 'geometry']]
        .rename(columns={'name': 'Country'})
    )
    gpkg_path = os.path.join(BASE_OUTPUT_DIR, 'country_summary.gpkg')
    gpd.GeoDataFrame(
        grouped.merge(countries_4326, on='Country', how='left'),
        crs='EPSG:4326',
    ).to_file(gpkg_path, driver='GPKG')
    print('  Saved: {}'.format(os.path.basename(gpkg_path)))

    _print_summary(agg_gdf, grouped)
    print('\nDone in {:.2f} minutes'.format((time.time() - t0) / 60))
