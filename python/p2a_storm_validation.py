'''
IKI Bangladesh (MIOASI): Storm verification plots

This script compares downscaled wind and pressure timeseries of 9 emsemble runs with
ERA5 and IBTrACS data.  It uses dask for parallelisation of storms.

Example slurm submission:

#!/bin/bash -l
#SBATCH --qos=normal
#SBATCH --mem=50G
#SBATCH --ntasks=14
#SBATCH --time=00-01:00:00
#SBATCH --export=NONE

python3 -u p2a_storm_validation.py

Note: Python 3 compatible only.

Author: ZM & HS
Created: 25/12/19
QA: HS
'''
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
import multiprocessing
import numpy as np
import dask
import dask.bag as db
import dataprocessing as dp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
from ascend import shape
from pandas.plotting import register_matplotlib_converters
from user_vars import ERA5DIR, EVENTS, HCNC


def process_storm(dfrow):
    """
    Process one storm (row)
    Args:
        dfrow: EVENTS dataframe row

    Returns: None

    """
    print(f'Processing {dfrow.NAME}... ')
    # loading the downsclaed data
    gust44 = dp.get_ds_storm_ts(HCNC, 'fg.T1Hmax', dfrow.NAME, RES, shpmask=val_shape)
    pres44 = dp.get_ds_storm_ts(HCNC, 'psl.T1Hmin', dfrow.NAME, RES, shpmask=val_shape)
    pres44 = pres44/100  # Convert to hPa
    wind44 = dp.get_wspd_ts(HCNC, dfrow.NAME, RES, shpmask=val_shape)
    # converting gust units m/sec to knots
    gust44 = dp.mpsec_to_knots(gust44)
    wind44 = dp.mpsec_to_knots(wind44)
    # Define time bounds for plotting
    starttime = gust44.index[0] - pd.Timedelta(days=DAYS)  # start time 3 days before start fo downscaled data
    endtime = gust44.index[-1]  # Get end time from last step of last run of downscaled
    # loading ibtracs data
    ibtracs = dp.get_ibtracs_ts(IBTRACSFILE, dfrow.IBID)[starttime:endtime]
    werror = 10
    perror = 3
    # loading ERA5 data
    era5 = dp.get_era5_ts(ERA5DIR, dfrow.NAME, shpmask=val_shape)[starttime:endtime]
    era5.PRES = era5.PRES/100  # Convert to hPa


    # Plotting
    gust_color_cycle = ['#edf8b1', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#253494', '#081d58',
                        '#000000']
    wind_color_cycle = ['#fff7f3', '#fde0dd', '#fcc5c0', '#fa9fb5', '#f768a1', '#dd3497', '#ae017e', '#7a0177',
                        '#49006a']
    psl_color_cycle = ['#ccece6', '#99d8c9', '#66c2a4', '#41ae76', '#238b45', '#006d2c', '#00441b', '#002d12',
                       '#000000']

    fig, axs = plt.subplots(3, 1, figsize=(7, 7.5), sharex='col')
    fig.suptitle(f'{dfrow.NAME}', fontsize=12, weight='bold')
    for ax in [axs[0], axs[1]]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    axs[2].spines['top'].set_visible(False)
    axs[2].spines['right'].set_visible(False)
    # Axis labels
    axs[0].set_ylabel('Speed (knots)')
    axs[1].set_ylabel('Speed (knots)')
    axs[2].set_ylabel('Pressure (hPa)')
    for ax in axs:
        ax.xaxis_date()
        ax.xaxis.set_major_locator(mdates.HourLocator(byhour=[0, 12]))
        ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=[3, 6, 9, 15, 18, 21]))
    # Axis labels
    axs[0].text(0.03, 0.9, 'Maximum 3-sec Gust Speed', color='darkgrey', fontsize=10, ha='left', weight='bold',
                transform=axs[0].transAxes)
    axs[1].text(0.03, 0.9, 'Maximum Wind Speed', color='darkgrey', fontsize=10, ha='left', weight='bold',
                transform=axs[1].transAxes)
    axs[2].text(0.03, 0.1, 'Minimum Mean \nSea-level Pressure', color='darkgrey', fontsize=10, ha='left', va='bottom', weight='bold',
                transform=axs[2].transAxes)
    # Gust speed plotting
    l5 = axs[0].plot(era5.index, era5.GUST, 'k-', lw=2, zorder=200)
    # axs[0].plot(era5.index, era5.WIND, 'r-', lw=1, zorder=200)
    
    for i in range(9):
        l3 = axs[0].plot(gust44.index, gust44[i], color=gust_color_cycle[i], lw=3, label='Gust Speed')
        axs[0].plot(gust44.index, gust44[i], color='white', lw=1)
    # Wind Speed Plotting
    axs[1].errorbar(ibtracs.index - pd.Timedelta(minutes=20), ibtracs.USA_WIND, yerr=werror, fmt='<', ms=4,
                    zorder=100, label="IBTrACS US", color='grey', mfc='darkgrey', lw=1)
    axs[1].errorbar(ibtracs.index + pd.Timedelta(minutes=20), ibtracs.NEWDELHI_WIND, yerr=werror, fmt='>', ms=4,
                    zorder=100, label="IBTrACS ND", color='grey', mfc='darkgrey', lw=1)
    axs[1].plot(era5.index, era5.WIND, 'k-', lw=2, zorder=200)
    for i in range(9):
        l4 = axs[1].plot(wind44.index, wind44[i], color=wind_color_cycle[i], lw=2, label='Wind Speed')
        axs[1].plot(wind44.index, wind44[i], color='white', lw=1)
    
    # MSLP plotting
    l5 = axs[2].plot(era5.index, era5.PRES, 'k-', lw=2, label='ERA5', zorder=200)
    # axs[1].text(endtime + pd.Timedelta(hours=1), era5.PRES[-1] / 100, 'ERA5', va='center', weight='bold')
    axs[2].errorbar(ibtracs.index - pd.Timedelta(minutes=20), ibtracs.USA_PRES, yerr=perror, fmt='<', ms=4,
                    zorder=100, label="IBTrACS US", color='grey', mfc='darkgrey', lw=1)
    axs[2].errorbar(ibtracs.index + pd.Timedelta(minutes=20), ibtracs.NEWDELHI_PRES, yerr=perror, fmt='>', ms=4,
                    zorder=100, label="IBTrACS ND", color='grey', mfc='darkgrey', lw=1)
    for i in range(9):
        l4 = axs[2].plot(pres44.index, pres44[i], lw=3, color=psl_color_cycle[i], label='_nolegend_')
        axs[2].plot(pres44.index, pres44[i], lw=1, color='white', label='_nolegend_')
    axs[2].legend(frameon=False, fontsize=8, loc=4)
    axs[2].xaxis.set_major_formatter(mdates.DateFormatter('%b-%d\n%HZ'))
    plt.tight_layout()
    # plt.show()
    fig.savefig(f'/project/ciid/projects/IKI/historical_catalogue/plots/compare{dfrow.NAME}+wind.png', dpi=600)
    plt.close()
    print(f'Done {dfrow.NAME}!')

# Pandas to matplotlib datetime conversion handling etc.
register_matplotlib_converters()

RES = '4p4'  # '4p4' or '1p5'
DAYS = 1  # Days pre downscale data to plot

IBTRACSFILE = '/project/ciid/projects/IKI/obs_datasets/global/IBTrACS/v04r00/ibtracs.NI.list.v04r00.csv'

# Load Validation area shapefile
val = shape.load_shp('/project/ciid/projects/IKI/ArcGIS/ValidationArea2.shp')
val_shape = val.unary_union()


# # Dask Parallel processing
# # Determine the number of processors visible...
# cpu_count = multiprocessing.cpu_count()
# # .. or as given by slurm allocation.
# if 'SLURM_NTASKS' in os.environ:
#     cpu_count = os.environ['SLURM_NTASKS']
# # Do not exceed the number of CPU's available, leaving 1 for the system.
# num_workers = cpu_count - 1
# print('Using {} workers from {} CPUs...'.format(num_workers, cpu_count))

with dask.config.set(num_workers=13):
    dask_bag = db.from_sequence([row for index, row in EVENTS.iterrows()]).map(process_storm)
    dask_bag.compute()
    # For debugging
    # dask_bag.compute(scheduler='single-threaded')
