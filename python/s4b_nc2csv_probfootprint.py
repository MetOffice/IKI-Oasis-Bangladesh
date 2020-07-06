'''
IKI Bangladesh (MIOASI): Data format conversion from netCDF to OASIS format (csv) for footprint files only!

This script specifically creates Oasis footprint files with non-1 probabilities.  All other Oasis files can be
re-combined (with some minor modificiations).

Note: Python 3 compatible only

Authors: HS
Created: 20/4/19
'''
import pandas as pd
import iris
import numpy as np
import pdb

NUMBINS = 110  # Numher of intensity bins dictating the 3rd dimension of the netCDF files

events = ['fp.exceedance2.fg.T1Hmax.AILA.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.AKASH.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.BOB01.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.BOB07.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.BULBUL.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.FANI.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.MORA.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.RASHMI.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.ROANU.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.SIDR.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.TC01B.4p4km.nc',
          'fp.exceedance2.fg.T1Hmax.VIYARU.4p4km.nc']

footprintdf = pd.DataFrame(columns=['EVENT_ID', 'AREAPERIL_ID', 'INTENSITY_BIN_INDEX', 'PROBABILITY'])
# Set up CSV file
outfile = 'hcprob_oasis_footprint.csv'
footprintdf.to_csv(outfile, mode='w', index=None, header=True)

for enum, event in enumerate(events):
    cube = iris.load_cube(event)
    # diffcubedata = np.abs(np.diff(cube.data.data))
    diffcubedata = cube.data.data
    for binid in range(NUMBINS):
        datasize = diffcubedata[:, :, binid].size
        eventIDs = np.ones(datasize) * (enum + 1)
        binsliceINDEXs = np.ones(datasize) * (binid + 1)
        
        footprintdf = pd.DataFrame({'EVENT_ID': eventIDs.astype('int'),
                                    'AREAPERIL_ID': np.array(range(datasize)) + 1,
                                    'INTENSITY_BIN_INDEX': binsliceINDEXs.astype('int'),
                                    'PROBABILITY': diffcubedata[:, :, binid].flatten()})
        # Remove lines with 0 probability
        footprintdf = footprintdf[footprintdf.PROBABILITY != 0]
        # Append to csv file
        footprintdf.to_csv(outfile, mode='a', index=None, header=False)

# Load and sort
df = pd.read_csv(outfile)
# Sort columns
df.sort_values(by=['EVENT_ID', 'AREAPERIL_ID', 'INTENSITY_BIN_INDEX'], inplace=True)
df.astype({'PROBABILITY': np.float64}).to_csv(outfile, mode='w', index=None, header=True)
df.to_csv(outfile, mode='w', index=None, header=True, float_format='%.7f')