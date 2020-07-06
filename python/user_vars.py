'''
IKI Bangladesh (MIOASI): User vars & constants

User variables and constants for use in the IKI Bangladesh project.

Note: Python 3 compatible only

Author: Hamish Steptoe (hamish.steptoe@metoffice.gov.uk)
Created: 7/3/19
QA:
'''
import pandas as pd

# DIRECTORIES
HCPP = '/project/ciid/projects/IKI/historical_catalogue/downscaled/pp'
SCRATCH = '/spice/scratch/hsteptoe/iki'
HCNC = '/project/ciid/projects/IKI/historical_catalogue/downscaled/netcdf'
ESDIR = '/project/ciid/projects/IKI/eventset'
ERA5DIR = '/project/ciid/projects/IKI/era5fps'

DOMAINS = {  # Name: [lon min, lon max, lat min, lat max]
    'BGD': [87.5, 93.0, 20.5, 27.5]  # TODO: Check this domain is up to date with the latest downscaling suite domain
}

# Define storms, and initial and final cycle points
# This is useful for globbing pp files
# Set event-name and cycle points to get correct pp files later
EVENTS = pd.read_csv('/net/home/h03/hsteptoe/lib/IKI/trunk/events.csv', header=4)


# Resolution options
RES = ['4p4', '1p5']

# Dict for translating long name to short name for filenaming purposes
NAMES = {'net_down_surface_sw_flux_corrected': 'rsnds',
         'wet_bulb_potential_temperature': 'wbpt',
         'turbulent_kinetic_energy': 'tke',
         'air_pressure_at_sea_level': 'psl',
         'air_temperature': 'tas',
         'geopotential_height': 'zg',
         'relative_humidity': 'hur',
         'stratiform_rainfall_amount': 'prlst',
         'stratiform_snowfall_amount': 'prlssn',
         'surface_downwelling_shortwave_flux_in_air': 'rsds',
         'wind_speed_of_gust': 'fg',
         'x_wind': 'ua',
         'y_wind': 'va'}

# TC Intensity Scale for North Indian basin
# see http://www.wmo.int/pages/prog/www/tcp/operational-plans.html
tci = [0, 9, 14, 17, 25, 33, 46, 62]
tcin = ['', 'Depression', 'Deep Depression', 'Cyclonic storm', 'Severe cyclonic storm',
        'Very severe cyclonic storm', 'Extremely severe cyclonic storm',
        'Super cyclonic storm']
