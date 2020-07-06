'''
IKI Bangladesh (MIOASI): S1b Tidy netCDF metadata

In some instances, it's useful to run this script independent of other
data processing scripts.

Author: HS
Created: 19/7/19
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
from user_vars import EVENTS, NAMES, RES, SCRATCH, HCNC

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
ATTRS = {'summary': 'Tropical cyclone footprint ensemble over Bangladesh',
         'title': 'Downscaled Tropical Cyclone footprint ensemble over Bangladesh',
         'date_created': time.strftime('%Y%m%dT%H:%M:%S'),
         'contact': 'enquiries@metoffice.gov.uk',
         'Conventions': 'CF-1.7',  # Name of the conventions followed by the dataset.
         'comment': 'Supported by the International Climate Initiative (IKI) and the Federal Ministry for the '
                    'Environment, Nature Conservation and Nuclear Safety, based on a decision of the Germany '
                    'Bundestag',
         'data_type': 'grid',
         'spatial_resolution': pargs.resolution.replace('p', '.') + 'km',
         'history': '(1.0) Initial release',
         'keywords': 'Bangladesh, ensemble, footprint, Met Office',
         'product_version': 'v1.0',
         'project': 'Oasis Platform for Climate and Catastrophe Risk Assessment – Asia',
         'references': '',  # References that describe the data or methods used to produce it
         'source': 'Copernicus Climate Change Service Information (C3S) ECMWF ERA5 / Met Office UM RA2T CON',
         # Method of production of the original data.
         'standard_name_vocabulary': 'NetCDF Climate and Forecast (CF) Standard Names version 51',
         'type': 'float',
         'licence': 'Creative Commons Attribution 4.0 International (CC BY 4.0)'
         }

dp.tidy_netcdf(HCNC, fglob='fpens*.nc', global_attributes=ATTRS)
