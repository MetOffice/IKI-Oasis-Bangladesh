'''
IKI Bangladesh (MIOASI): Calculate min/max dataset differences

Calculate min/max values of ERA5, IBTrACS and RA2 datasets for gust or MSLP.
Uses dask to parallelise the process. Run on SPICE.

#!/bin/bash -l
#SBATCH --qos=normal
#SBATCH --mem=20G
#SBATCH --ntasks=14
#SBATCH --time=00-00:10:00
#SBATCH --export=NONE

python3 -u s7b_validation_difference2.py

Python 3 compatible only.

Author: HS
Created: 26/5/20
'''
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
import numpy as np
import dask
import dask.bag as db
import dataprocessing as dp
import matplotlib as mpl

mpl.use('Agg')
import pandas as pd
from ascend import shape
from pandas.plotting import register_matplotlib_converters
from user_vars import ERA5DIR, EVENTS, HCNC
import pdb


def _make_dataframe(model, storm, var, inten, time):
    df = pd.DataFrame({'MODEL': model,
                       'STORM': storm,
                       'VAR': var,
                       'INTEN': inten,
                       'TIME': time})
    return df


def matrixmin(a, b):
    if len(b) == 0:
        return np.array([np.nan] * 9)
    # elif len(b) == 1:
    #     return a
    else:
        # Return the times of b that have the smallest difference with times of a
        xs, ys = np.meshgrid(a, b)
        return b[np.abs(xs - ys).argmin(axis=0)]


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
    pres44 = dp.get_ds_storm_ts(HCNC, 'psl.T1Hmin', dfrow.NAME, RES, shpmask=val_shape) / 100.  # Convert to hPa
    wind44 = dp.get_wspd_ts(HCNC, dfrow.NAME, RES, shpmask=val_shape)
    # converting gust units m/sec to knots
    gust44 = dp.mpsec_to_knots(gust44)
    wind44 = dp.mpsec_to_knots(wind44)
    # Define time bounds
    starttime = gust44.index[0]  # start time of downscaled data
    endtime = gust44.index[-1]  # Get end time from last step of last run of downscaled
    
    # loading ibtracs data
    ibtracs = dp.get_ibtracs_ts(IBTRACSFILE, dfrow.IBID)[starttime:endtime]
    
    # loading ERA5 data
    era5 = dp.get_era5_ts(ERA5DIR, dfrow.NAME, shpmask=val_shape)[starttime:endtime]
    era5.PRES = era5.PRES / 100  # Convert to hPa

    # Find absolute values: ERA5
    era5dfg = _make_dataframe('ERA5', dfrow.NAME, 'gust', [era5.GUST.max()]*9, era5.GUST.idxmax())
    era5dfw = _make_dataframe('ERA5', dfrow.NAME, 'wind', [era5.WIND.max()]*9, era5.WIND.idxmax())
    era5dfp = _make_dataframe('ERA5', dfrow.NAME, 'mslp', [era5.PRES.min()]*9, era5.PRES.idxmin())
    
    # RA2
    ra2dfg = _make_dataframe('RA2', dfrow.NAME, 'gust', gust44.max(), gust44.idxmax())
    ra2dfw = _make_dataframe('RA2', dfrow.NAME, 'wind', wind44.max(), wind44.idxmax())
    ra2dfp = _make_dataframe('RA2', dfrow.NAME, 'mslp', pres44.min(), pres44.idxmin())
    
    # Account for possible join max/min times in IBTrACS data
    # Resample to match RA2 hourly output
    ibndw = ibtracs[ibtracs.NEWDELHI_WIND == ibtracs.NEWDELHI_WIND.max()].resample('1H').pad().index.values
    ibusw = ibtracs[ibtracs.USA_WIND == ibtracs.USA_WIND.max()].resample('1H').pad().index.values
    
    ibndp = ibtracs[ibtracs.NEWDELHI_PRES == ibtracs.NEWDELHI_PRES.min()].resample('1H').pad().index.values
    ibusp = ibtracs[ibtracs.USA_PRES == ibtracs.USA_PRES.min()].resample('1H').pad().index.values
    
    # New Delhi absolute values
    nddfw = _make_dataframe('IBND', 
                            dfrow.NAME, 
                            'wind', 
                            ibtracs.NEWDELHI_WIND.max(),
                            matrixmin(wind44.idxmax().values, ibndw))

    nddfp = _make_dataframe('IBND', 
                            dfrow.NAME, 
                            'mslp', 
                            ibtracs.NEWDELHI_PRES.min(),
                            matrixmin(pres44.idxmax().values, ibndp))
    
    # US absolute values
    usdfw = _make_dataframe('IBUS', 
                            dfrow.NAME, 
                            'wind', 
                            ibtracs.USA_WIND.max(),
                            matrixmin(wind44.idxmax().values, ibusw))
    
    usdfp = _make_dataframe('IBUS', 
                            dfrow.NAME, 
                            'mslp', 
                            ibtracs.USA_PRES.min(),
                            matrixmin(pres44.idxmin().values, ibusp))
    
    # Return differnce dataframes
    print(f'Done {dfrow.NAME}!')
    return pd.concat([usdfw, usdfp, nddfw, nddfp, era5dfg, era5dfw, era5dfp, ra2dfg, ra2dfw, ra2dfp], ignore_index=True)


# Pandas to matplotlib datetime conversion handling etc.
register_matplotlib_converters()

RES = '4p4'  # '4p4' or '1p5'
DAYS = 1  # Days pre downscale data to plot

IBTRACSFILE = 'sup/ibtracs.NI.list.v04r00.csv'

# Load Validation area shapefile
val = shape.load_shp('sup/ValidationArea2.shp')
val_shape = val.unary_union()

# Start csv file for recording differences
DIFFCSV = 'validation_diff2.csv'

# # Dask Parallel processing
# # Determine the number of processors visible...
# cpu_count = multiprocessing.cpu_count()
# # .. or as given by slurm allocation.
# if 'SLURM_NTASKS' in os.environ:
#     cpu_count = os.environ['SLURM_NTASKS']
# # Do not exceed the number of CPU's available, leaving 1 for the system.
# num_workers = cpu_count - 1
# print('Using {} workers from {} CPUs...'.format(num_workers, cpu_count))

# EVENTS = EVENTS[EVENTS.NAME != "FANI"]

with dask.config.set(num_workers=13):
    dask_bag = db.from_sequence([row for index, row in EVENTS.iterrows()]).map(process_storm)
    # dfs = dask_bag.compute(scheduler='single-threaded')
    dfs = dask_bag.compute()

diffdf = pd.concat(dfs, ignore_index=True)
diffdf.sort_values(by=['VAR', 'STORM', 'MODEL'], inplace=True)
# Write to CSV
diffdf.to_csv(DIFFCSV, mode='w', index=None, header=True)
