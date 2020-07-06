'''
IKI Bangladesh (MIOASI): Data format conversion from netCDF to OASIS format (csv)

Note: Python 3 compatible only

Authors: BB, updated by ZM
Created: 12/12/19
QA: HS
'''
import glob
import iris
import numpy as np
import os
import pandas as pd
from user_vars import DOMAINS, HCNC

# Original netCDF file directory
INPUTDIR = HCNC
# File pattern
FTYPE = 'fpens'  # Identifier for file type
VAR = 'fg.T1Hmax'  # Variable and time meaning
RES = '4p4'  # Data resolution (4p4 or 1p5)
FILEPAT = f'{FTYPE}*{VAR}*[!grand].{RES}*.nc'
# Dimension that defines the ensemble number
# ENSD = None   # ...if there are no defined ensemble members
ENSD = 'forecast_reference_time'
# Output directory
OUTDIR = f'/oasis_csv/{RES}'
OUTPREFIX = 'hc'  # Output file prefix, hc = historical catalogue, es = event set
# Global bin vars
RANGEMIN = 0  # Minimum expected range of input files
RANGEMAX = 105  # Maximum expected range of input files
BINWIDTH = 0.1  # Resolving bin width for data


def write_df(event, index, aaaa, bbbb, header, id, cccc=None):
    """

     write out an # x, y, z file
     output could be event_id, index, aaaa,      bbbb,     for all grid points:
     or              event_id, index, aaaa,      bbbb,  cccc,    for all grid points:
     i.e.:
     output could be event_id, index, longitude, latitude, for all grid points:

    """
    # print event
    no_nan = ~np.isnan(aaaa)
    no_nan_flat = no_nan.flatten(order='F')
    # print(no_nan, no_nan.shape, no_nan_flat.shape)
    event = event[no_nan_flat]
    index = index[no_nan_flat]
    # print(np.shape(index), index, ' = index! ')
    aaaa = aaaa[no_nan]
    # bbbb may be a complete 2-D field with nans, such as latitudes or a 1-D field without nans such as confidence
    try:
        bbbb = bbbb[no_nan_flat]
    except IndexError:
        bbbb = bbbb[no_nan]
    
    if cccc is None:
        dat = np.concatenate([event.flatten(order='F')[..., np.newaxis],
                              index.flatten(order='F')[..., np.newaxis],
                              aaaa.flatten(order='F')[..., np.newaxis],
                              bbbb.flatten(order='F')[..., np.newaxis]],
                             axis=1)
    else:
        cccc = cccc[no_nan_flat]
        dat = np.concatenate([event.flatten(order='F')[..., np.newaxis],
                              index.flatten(order='F')[..., np.newaxis],
                              aaaa.flatten(order='F')[..., np.newaxis],
                              bbbb.flatten(order='F')[..., np.newaxis],
                              cccc.flatten(order='F')[..., np.newaxis]],
                             axis=1)
    
    df = pd.DataFrame(data=dat, columns=header)
    
    if id == 'sfp':
        df = df.astype({header[0]: int, header[1]: int, header[2]: int, header[3]: float})
    elif id == 'latlon':
        df = df.astype({header[0]: int, header[1]: int, header[2]: float, header[3]: float})
    elif id == 'bin':
        df = df.astype({header[0]: int, header[1]: float, header[2]: float, header[3]: float, header[3]: float})
    
    return df


def write_latlon(event_number, cube, index, oasis_file):
    """

     write out one  oasis compatible csv files
     write out  one file for locations

     output 1. should be index, areaperilid, longitude, latitude, for all grid points:

     read  the grid cartesian coordinates
    """
    # Check we get correct cube dimension for function
    assert cube.ndim == 2, 'Cube should have 2 dimension coordinates only ie. be 2-D '
    
    lon = cube.coord('longitude').points
    lat = cube.coord('latitude').points
    
    lons, lats = np.meshgrid(lon, lat)
    
    lats = lats[~np.isnan(cube.data)]
    lons = lons[~np.isnan(cube.data)]
    
    header = ['EVENT_ID', 'AREAPERIL_ID', 'LONGITUDE', 'LATITUDE']
    id = 'latlon'
    
    df = write_df(event_number, index, lons, lats, header, id)
    write_csvfile(df, oasis_file)
    
    return None


def write_period_id(df):
    occyear = df['OCC_YEAR']
    periodid = (occyear - min(occyear)) + 1
    df['PERIOD_NO'] = periodid
    
    return df


def write_period(eventid, cube, df):
    """

     returns the dataframe for event occurance time

     output  should be  event_id, period_no, occ_year, occ_month, occ_day

     Name 	Type 	Bytes 	Description 	                                Example
     event_id 	int 	4 	The occurrence event_id (for every ensemble)    45567
     period_no 	int 	4 	A numbered period in which the event occurs 	56876
     occ_year 	int 	4 	the year number of the event occurrence 	    1999
     occ_month 	int 	4 	the month of the event occurrence               5
     occ_day 	int 	4 	the day of the event occurrence 		        16

    """
    if ENSD is not None:
        time = cube.coord("forecast_reference_time")
        dates = time.units.num2date(time.points)
        
        yr = dates[-1].year
        mm = dates[-1].month
        dd = dates[-1].day
        pr = 0
        
        df = df.append({'EVENT_ID': eventid, 'PERIOD_NO': pr, 'OCC_YEAR': yr, 'OCC_MONTH': mm, 'OCC_DAY': dd},
                       ignore_index=True)
    else:  # TODO: Check this is the best way to handle the case where ENSD = None
        print('No period coordinate found, setting period to 1/1/1!')
        df = df.append({'EVENT_ID': count, 'PERIOD_NO': 1, 'OCC_YEAR': 1, 'OCC_MONTH': 1, 'OCC_DAY': 1},
                       ignore_index=True)
    
    return df


def write_values(event_number, cube, index, oasis_file):
    """

     write out one  oasis compatible csv file
     write out one file  for wind gusts histogram bin index an # x, y, z file

     output  should be index, grid point histogram bin index value, confidence   for all grid points:

    """
    # Check we get correct cube dimension for function
    assert cube.ndim == 2, 'Cube should have 2 dimension coordinates only ie. be 2-D '
    
    # Confidence set to 1 for all grid point histogram bin index values
    confidence = np.ones(cube.data.size)
    
    # Convert cube.data to an index in a histogram
    # Histogram has BINWIDTH bin widths (eg. 0.1 m/s), eg. 900 bins for a value range between 0 and 90 m/s
    cube = np.ma.masked_invalid(cube.data)
    if np.any(cube.data > RANGEMAX):
        raise ValueError(f'Cube data range [{np.max(cube.data)}] exceeds max specified data range [{RANGEMAX}]!\n'
                         f'Adjust RANGEMAX and re-run...')
    gustbins = np.floor(cube.data / BINWIDTH)
    
    header = ['EVENT_ID', 'AREAPERIL_ID', 'INTENSITY_BIN_INDEX', 'PROBABILITY']
    id = 'sfp'
    
    df = write_df(event_number, index, gustbins, confidence, header, id)
    write_csvfile(df, oasis_file)


def write_occurrence(eventid, ensemblecount, filename, df):
    """
     returns the dataframe for mapping event id and event occurnace id
    """
    description = filename.split(".")[3] + ' ensemble ' + str(ensemblecount)
    df = df.append({'EVENT_ID': eventid, 'EVENT_OCCURANCE_ID': eventid, 'DESCRIPTION': description},
                   ignore_index=True)
    
    return df


def write_eventid(eventid, df):
    """
     output  should be  event_id

     Name 	    Type 	Bytes 	Description 	                Example
     event_id 	int 	4 	   The occurrence event_id 	    45567

    """
    df = df.append({'EVENT_ID': eventid}, ignore_index=True)
    
    return df


def write_csvfile(df, oasis_out):
    """
     write out oasis compatible csv files
    """
    if os.path.isfile(oasis_out):
        df.to_csv(oasis_out, mode='a', index=None, header=False)
    else:
        df.to_csv(oasis_out, mode='w', index=None, header=True)


def write_bindict(oasis_file):
    """
     write out one  oasis compatible csv files
     write out one file  for bin dictionary       # x, y, z file

     resolving range [ 0. - 90. ] in 0.1 m/s bins

     output  should be  bin index, bin_from, bin_to, interpolation, interval_type
    """
    #   needs a bin dictionary a.k.a.
    #   write once for all storm foot prints
    #   BIN_INDEX, BIN_FROM, BIN_TO, INTERPOLATION, INTERVAL_TYPE
    bin_index = np.arange((RANGEMAX - RANGEMIN) / BINWIDTH)
    bin_index = bin_index + 1
    bin_from = (bin_index - 1) * BINWIDTH
    bin_to = bin_from + BINWIDTH
    delta = 0.5 * BINWIDTH
    bin_mid = bin_from + delta
    confidence = np.ones(int(np.max(bin_index)))
    interval_type = confidence * 1202  # 1202 some code for OASIS compliance
    
    header = ['BIN_INDEX', 'BIN_FROM', 'BIN_TO', 'INTERPOLATION', 'INTERVAL_TYPE']
    
    id = 'bin'
    df = write_df(bin_index, bin_from, bin_to, bin_mid, header, id,
                  cccc=interval_type)
    write_csvfile(df, oasis_file)


def main():
    """
        standalone test program
        Input:  None
    """
    # Check OUTDIR
    if not os.path.exists(OUTDIR):
        print(f'Making directory: {OUTDIR}...')
        os.makedirs(OUTDIR)
    else:
        # Check if file exist.  Delete if they do as write functions make use of append.
        print(f'Checking for old files in {OUTDIR}...')
        for f in glob.glob(f'{OUTDIR}/*.csv'):
            os.remove(f)
            
    # Get original netCDF files, excluding any files with 'grand' in their name
    FILEPAT = f'{FTYPE}*{VAR}*[!grand].{RES}*.nc'
    all_files = sorted(glob.glob(f'{INPUTDIR}/{FILEPAT}'))
    
    eventoccur_df = pd.DataFrame(columns=['EVENT_ID', 'EVENT_OCCURANCE_ID', 'DESCRIPTION'])
    period_df = pd.DataFrame(columns=['EVENT_ID', 'PERIOD_NO', 'OCC_YEAR', 'OCC_MONTH', 'OCC_DAY'])
    eventid_df = pd.DataFrame(columns=['EVENT_ID'])
    
    # Setup unique EVENTID counter
    EVENTID = 1
    for filenum, filename in enumerate(all_files):
        # Load each cube
        print(f'Processing {filename}...')
        cube = iris.load_cube(filename)
        # cube = cube.intersection(longitude=(DOMAINS["BGD"][0], DOMAINS["BGD"][1]),
        #                          latitude=(DOMAINS["BGD"][2], DOMAINS["BGD"][3]))
        
        # Check for ensemble members based on ENSD dimension
        if ENSD is None:
            assert cube.ndim == 2, f'Cube should have 2 dimensions, but it has {cube.ndim}'
            n_ensembles = 1  # TODO: Would this work for the eventset data
        elif ENSD in [x.name() for x in cube.coords(dim_coords=True)]:
            n_ensembles = cube.coord(ENSD).shape[0]
        else:
            raise ValueError(f'Dimension {ENSD} does not exist in cube!')
            
        # df = pd.DataFrame(columns=['EVENT_ID', 'EVENT_OCCURANCE_ID', 'DESCRIPTION'])
        for ensnum in range(n_ensembles):
            print(f'Writing Eventid: {EVENTID} {filename} ensemble {ensnum+1}')

            if n_ensembles > 1:
                ecube = cube[ensnum, :, :]
            else:
                ecube = cube
 
            n_elem = ecube.data.size
            
            # catalogue_event =  event occurance id containing all emsembles
            catalogue_event = np.ones(n_elem) * EVENTID
            
            # index: used as areaperil_id
            index = np.array(range(n_elem)) + 1
            
            period_df = write_period(EVENTID, cube, period_df)
            eventoccur_df = write_occurrence(EVENTID, (ensnum+1), filename, eventoccur_df)
            write_latlon(catalogue_event, ecube, index, f'{OUTDIR}/{OUTPREFIX}_oasis_locations.csv')
            write_values(catalogue_event, ecube, index, f'{OUTDIR}/{OUTPREFIX}_oasis_footprint.csv')
            eventid_df = write_eventid(EVENTID, eventid_df)
            EVENTID += 1
    
    period_df = write_period_id(period_df)
    write_bindict(f'{OUTDIR}/{OUTPREFIX}_oasis_damage_bin_dict.csv')
    write_csvfile(eventid_df, f'{OUTDIR}/{OUTPREFIX}_oasis_events.csv')
    write_csvfile(period_df, f'{OUTDIR}/{OUTPREFIX}_oasis_occurrence.csv')
    write_csvfile(eventoccur_df, f'{OUTDIR}/{OUTPREFIX}_oasis_occurrenceid.csv')
    # Update dir and file permissions
    os.chmod(OUTDIR, 0o775)
    for f in glob.glob(f'{OUTDIR}/*.csv'):
        os.chmod(f, 0o775)
    


#
# MODULE OR STAND-ALONE
# _________________________________________________________________________

if __name__ == "__main__":
    main()
