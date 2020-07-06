'''
IKI Bangladesh (MIOASI): S1 Make 'basic' netcdf files from pp files

Make netcdf files from the raw pp file output. 'Basic' in the sense that
no further post-processing is performed.  This script mainly manipulates
netcdf metadata and splits the pp files into separate nc files.

Python 3 only.

Usage:
    Recommended to run in SPICE, e.g.:
    
    #!/bin/bash -l
    #SBATCH --qos=normal
    #SBATCH --mem=10G
    #SBATCH --ntasks=2
    #SBATCH --time=00-04:00:00
    #SBATCH --export=NONE
    module load scitools
    python3 -u s1_make_basic_netcdf.py -r 4p4 -e Fani -t

Author: HS
Created: 19/7/19
QA: NS 5/8/2019
'''

import argparse
import datetime as dt
import glob
import iris
import numpy as np
import os
import sys
import time
from cf_units import Unit

from ciid_tools.rim_remove import rim_remove
import dataprocessing as dp
from user_vars import EVENTS, NAMES, RES, SCRATCH

# Parse other processing options from the command line
parser = argparse.ArgumentParser(description='Specify file options')
parser.add_argument('-r', '--resolution', choices=RES, default='4p4',
                    help='Specify pp file resolution')
parser.add_argument('-e', '--event', choices=list(EVENTS.keys()), default='Sidr',
                    help='Specify event (storm name)')
parser.add_argument('-t', '--tidync', action='store_true',
                    help='Optional, run tidy netcdf function to tidy additional netcdf metadata')
pargs = parser.parse_args()

# Define constants
DATADIR = SCRATCH
OUTDIR = SCRATCH
RIM_WIDTH = 13

# Set custom global attributes specific to this dataset
ATTRS = {'summary': 'Tropical cyclone data over Bangladesh downscaled using Met Office RA2T_CON initiated from ERA5',
         'title': 'Downscaled Tropical Cyclone data over Bangladesh',
         'date_created': time.strftime('%Y%m%dT%H:%M:%S'),
         'contact': 'enquiries@metoffice.gov.uk',
         'Conventions': 'CF-1.7',  # Name of the conventions followed by the dataset.
         'comment': 'Supported by the International Climate Initiative (IKI) and the Federal Ministry for the '
                    'Environment, Nature Conservation and Nuclear Safety, based on a decision of the Germany '
                    'Bundestag',
         'data_type': 'grid',
         'spatial_resolution': pargs.resolution.replace('p', '.') + 'km',
         'history': '(1.0) Initial release',
         'keywords': 'Bangladesh, dynamical downscaling, RA2T, Met Office',
         'product_version': 'v1.0',
         'project': 'Oasis Platform for Climate and Catastrophe Risk Assessment – Asia',
         'references': '',  # References that describe the data or methods used to produce it
         'source': 'Copernicus Climate Change Service Information (C3S) ECMWF ERA5 / Met Office UM RA2T CON',
         # Method of production of the original data.
         'standard_name_vocabulary': 'NetCDF Climate and Forecast (CF) Standard Names version 51',
         'type': 'float',
         'licence': 'Creative Commons Attribution 4.0 International (CC BY 4.0)'
         }

# Get files
print(f"Making {pargs.event} {pargs.resolution}...")
# Glob all files
files = glob.glob(f"{DATADIR}/*{pargs.resolution}*.pp")
# Python can't do brace expansion, so then filter list to date range
# Make list of filters ie. start and end date to pass to range
flt = list(range(int(EVENTS[pargs.event]['ic'][:8]), int(EVENTS[pargs.event]['fc'][:8]) + 1))
# Convert back to strings
flt = [str(a) for a in flt]
# Do filter
files = [k for k in files for f in flt if f in k]
# Load into iris
cubes = iris.load(files)
# Amend cube metadata where needed
for c in cubes:
    if c.name() == 'm01s01i202':
        c.rename('net_down_surface_sw_flux_corrected')
        c.units = Unit('W m-2')
        if not c.cell_methods:
            c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
    elif c.name() == 'm01s03i473':
        # Optional - TKE is large c. 100GB on disk
        # c.rename('turbulent_kinetic_energy')
        # c.units = Unit('J kg-1')
        # if not c.cell_methods:
        #     c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
        continue
    elif c.name() == 'wet_bulb_potential_temperature' and not c.cell_methods:
        c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
    elif c.name() == 'air_pressure_at_sea_level' and not c.cell_methods:
        c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
    elif c.name() == 'air_temperature':
        if not c.cell_methods:
            c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
        if len(c.coord('forecast_period').points) == 2:
            if c.cell_methods[0].method == 'maximum':
                c.cell_methods = (iris.coords.CellMethod(method='maximum', coords='time'),)
            else:
                c.cell_methods = (iris.coords.CellMethod(method='minimum', coords='time'),)
    elif c.name() == 'geopotential_height':
        c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
    elif c.name() == 'relative_humidity' and not c.cell_methods:
        if len(c.coord('forecast_period').points) == 16:
            c.cell_methods = (iris.coords.CellMethod(method='point', coords='time', intervals='3 hour'),)
        else:
            c.cell_methods = (iris.coords.CellMethod(method='point', coords='time', intervals='1 hour'),)
    elif c.name() == 'surface_downwelling_shortwave_flux_in_air' and not c.cell_methods:
        c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
    elif c.name() == 'wind_speed_of_gust' and not c.cell_methods:
        c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
    elif c.name() == 'x_wind':
        if len(c.coord('forecast_period').points) == 16:
            c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
        else:
            c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
    elif c.name() == 'y_wind':
        if len(c.coord('forecast_period').points) == 16:
            c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
        else:
            c.cell_methods = (iris.coords.CellMethod(method='point', coords='time'),)
    # Some exceptions...
    elif c.name() in list(NAMES.keys()):
        # Some cubes don't need any modification
        # If we expect to see them, and they exist in NAMES,
        # carry on processing...
        pass
    else:
        # Otherwise it must be a name we aren't expecting
        # Mention in stdout but don't raise an error
        print(f'Cube name {c.name()} not recognised! Moving on...')
        # Move onto next c in for loop
        continue

    # Remove rim
    c = rim_remove(c, RIM_WIDTH)

    # Make time profile string
    if not c.coord('forecast_period').has_bounds():
        c.coord('forecast_period').guess_bounds()
    timediff = int(np.diff(c.coord('forecast_period').bounds[0]))
    tpro = 'T{}H{}'.format(timediff, c.cell_methods[0].method.replace('maximum', 'max').replace('minimum', 'min'))
    # Start date, IC point + 24 hours (we throw away first 24 hours as spin-up)
    sd = dt.datetime.strptime(EVENTS[pargs.event]['ic'], '%Y%m%dT%H%MZ') + dt.timedelta(hours=24)
    # End date, FC point + 72 hours
    ed = dt.datetime.strptime(EVENTS[pargs.event]['fc'], '%Y%m%dT%H%MZ') + dt.timedelta(hours=72)
    # Make save name
    outname = '{}.{}.UMRA2T.{}_{}.{}.{}km.nc'.format(NAMES[c.name()], tpro,
                                                     dt.datetime.strftime(sd, '%Y%m%d'),
                                                     dt.datetime.strftime(ed, '%Y%m%d'),
                                                     pargs.event.upper(), pargs.resolution)
    outfile = os.path.join(OUTDIR, outname)
    print(f'Saving {outfile}...')
    iris.save(c, outfile)
    
# Do tidy, read and overwrite files in OUTDIR
if pargs.tidync:
    dp.tidy_netcdf(OUTDIR, fglob='*UMRA2T*.nc', global_attributes=ATTRS)
