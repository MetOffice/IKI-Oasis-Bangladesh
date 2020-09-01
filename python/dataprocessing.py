'''
IKI Bangladesh (MIOASI): Data processing functions

Data processing functions relating to the IKI Bangladesh project.

Note: Python 3 compatible only

Author: HS
Created: 7/3/19
'''
import datetime
import math
import iris
import pandas as pd
import re
import time
import numpy as np
import netCDF4 as nc
import iris.coords as icoords
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
import glob
from six import integer_types
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import cascaded_union
from scipy.interpolate import interp1d
import pdb

def storm_name_callback(cube, field, filename):
    storm_name = re.split('\d{8}_\d{8}\.([^\.]+)', filename)[1]
    cube.attributes['Storm Name'] = storm_name.title()

def remove_forecast_coords_callback(cube, field, filename):
    try:
        cube.remove_coord('forecast_period')
        cube.remove_coord('forecast_reference_time')
    except:
        pass

def remove_forecast_period_callback(cube, field, filename):
    try:
        cube.remove_coord('forecast_period')
    except:
        pass

def lagged_ensemble_callback(cube, field, filename):
    # Add our own realization coordinate if it doesn't already exist.
    if not cube.coords('realization'):
        realization = np.int32(filename[-6:-3])
        ensemble_coord = icoords.AuxCoord(realization, standard_name='realization')
        cube.add_aux_coord(ensemble_coord)

def tidy_netcdf(dir, fglob=None, global_attributes=None, var_attributes=None):
    """Tidy and amend attributes of netCDF files

    Amend and add global attributes. Variable attributes can be set here also.

    Most of what is in 'attributes' part is required by the CF-1.7 convention
    but other stuff can be added if necessary.

    Arguments:
        dir (string): Location of netcdf files
        fglob (string): String for glob file names. Provides basic filtering functionality.
        global_attributes (dict): optonal dict of attributes to add or over-write default_attributes
        var_attributes (dict): optonal dict of variable specific attributes to
                               add or over-write exisitng variable attributes.

    Notes:
        All netcdf files within dir are given the same global attributes.
        The files seek to conform to CF-1.7 conventions
        http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7

    Returns:
        None. Netcdf files are modified in place (specified by dir).
    """
    if fglob is None:
        files = glob.glob(os.path.join(dir, '*.nc'))
    else:
        files = glob.glob(os.path.join(dir, fglob))

    for f in files:
        print("Processing {}".format(f))

        ncd = nc.Dataset(f, 'r+')
        # Set global attributes
        default_attributes = {'contact': 'enquiries@metoffice.gov.uk',
                              'Conventions': 'CF-1.7',  # Name of the conventions followed by the dataset.
                              'comment': '',  # Miscellaneous information about the data or methods used to produce it.
                              'data_type': 'grid',
                              'date_created': time.strftime('%Y%M%dT%H:%M:%S'),
                              'geospatial_lat_max': str(ncd.get_variables_by_attributes(axis='Y')[0][:] \
                                                    .max().astype('f')),
                              'geospatial_lat_min': str(ncd.get_variables_by_attributes(axis='Y')[0][:] \
                                                    .min().astype('f')),
                              'geospatial_lat_resolution': '{0:.2f}'.format(
                                                    ncd.get_variables_by_attributes(axis='Y')[0][1].astype('f') \
                                                    - ncd.get_variables_by_attributes(axis='Y')[0][0].astype('f')),
                              'geospatial_lat_units': 'degrees_north',
                              'geospatial_lon_max': str(ncd.get_variables_by_attributes(axis='X')[0][:] \
                                                    .max().astype('f')),
                              'geospatial_lon_min': str(ncd.get_variables_by_attributes(axis='X')[0][:] \
                                                    .min().astype('f')),
                              'geospatial_lon_resolution': '{0:.2f}'.format(
                                                    ncd.get_variables_by_attributes(axis='X')[0][1].astype('f') \
                                                    - ncd.get_variables_by_attributes(axis='X')[0][0].astype('f')),
                              'geospatial_lon_units': 'degrees_east',
                              'history': '(1.0) Initial release',
                              'id': os.path.basename(f),
                              'spatial_resolution': f[-8:-3].replace('p', '.'),
                              'institution': 'Met Office, UK',  # Where the data was produced
                              'keywords': '',
                              'product_version': '',
                              'project': '',
                              'references': '',  # References that describe the data or methods used to produce it
                              'source': '',  # Method of production of the original data.
                              'standard_name_vocabulary': 'NetCDF Climate and Forecast (CF) Standard Names version 68',
                              'summary': '',
                              'title': '',  # Short description of the file contents.
                              'type': ''}

        # If attributes is set, update default_attributes with new attributes
        if global_attributes is not None:
            for key in global_attributes.keys():
                default_attributes[key] = global_attributes[key]
        # Add attributes
        ncd.setncatts(default_attributes)

        # Set attributes specific to variables.
        # In some cases, the files contain multiple variables, but we only want to modify a specific one
        # Find intersection of netcdf variables and the ones we want to modify
        # This should catch instances of multiple vars in one netcdf file
        if var_attributes is not None:
            common_vars = set(ncd.variables.keys()) & set(var_attributes.keys())
            for att in list(common_vars):
                long_name = ncd.variables[att].long_name
                ncd.variables[att].setncatts(var_attributes[att][long_name])

        # Add latitude_longitude crs attributes
        coord_attributes = {'proj4': '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'}
        ncd.variables['latitude_longitude'].setncatts(coord_attributes)

        ncd.close()
        os.chmod(f, 0o774)


def calc_wspd(u, v):
    """
    function calculates the wind speed from wind components
    """
    return np.sqrt(u ** 2 + v ** 2)

def ws_units_func(u_cube, v_cube):
    if u_cube.units != getattr(v_cube, 'units', u_cube.units):
        raise ValueError("units do not match")
    return u_cube.units

def get_wspd_ts(path, storm, res, shpmask):
    """
    Extracts the U and V component and returns the wind speed timeseries of storm_dict

    Arguments:
        path (str): Path containing data to load
        storm (str): Name of storm
        res (str): Resolution of data

    Returns:
        Pandas dataframe with time index
    """
    ufile = f'{path}/ua.T1Hpoint.UMRA2T.*.{storm}.{res}km.nc'
    vfile = f'{path}/va.T1Hpoint.UMRA2T.*.{storm}.{res}km.nc'

    ucube = iris.load_cube(ufile, 'x_wind')
    vcube = iris.load_cube(vfile, 'y_wind')
    ucube = ucube.intersection(longitude=(75, 100), latitude=(10, 25))
    vcube = vcube.intersection(longitude=(75, 100), latitude=(10, 25))

    ws_ifunc = iris.analysis.maths.IFunc(calc_wspd, ws_units_func)
    ws_cube = ws_ifunc(ucube, vcube, new_name='wind speed')

    try:
        mwspd = shpmask.mask_cube(ws_cube)
    except:
        print("Can't mask with shape! Masked over lon-lat box instead...")
        mwspd = ws_cube
    
    cubedata = []
    timedata = []
    for subcube in mwspd.slices_over('forecast_reference_time'):
        # extracting the time
        tcoord = subcube.coord('time')
        units = tcoord.units
        tdata = [units.num2date(point) for point in tcoord.points]
        cube = subcube.collapsed(['latitude', 'longitude'], iris.analysis.MAX)

        cubedata.append(cube.data.filled())
        timedata.append(tdata)

    # Convert to Pandas Dataframe with unified time index
    s = list()
    [s.append(pd.Series(data=cubedata[i], index=timedata[i])) for i in range(np.shape(timedata)[0])]
    return pd.DataFrame(s).T


def get_ibtracs_ts(ifile, ibid):
    """
    Reads the ibtracs v04 csv file and extract wind and pressure time series of given storm

    Arguments:
        ifile (str): Path to IBTrACS file
        ibid (str): IBTrACS ID of storm

    Returns:
        Pandas dataframe
    """
    data = pd.read_csv(ifile,
                       index_col='ISO_TIME',
                       usecols=['SID', 'NAME', 'ISO_TIME', 'USA_WIND', 'USA_PRES', 'NEWDELHI_WIND', 'NEWDELHI_PRES','USA_LAT','USA_LON','NEWDELHI_LAT','NEWDELHI_LON'],
                       skiprows=[1], parse_dates=['ISO_TIME'], na_values=' ')

    return data[data.SID == ibid]


def get_ds_storm_ts(indir, var, storm, res, shpmask):
    """
    Load the downscaled historical_catalogue data (1.5 or 4.4 resolution) and return max/min pressure/wind.

    Arguments:
        indir (str): Path to data to load
        var (str): Full variable name, including time meaning eg. fg.T1Hmax
        storm (str): Storm name
        res (str): Data resolution

    Returns:
        Pandas dataframe with time index containing gust AND psl data
    """

    fpath = f'{indir}/{var}.*.{storm}.{res}km.nc'
    cube = iris.load_cube(fpath)

    cube = cube.intersection(longitude=(75, 100), latitude=(10, 25))
    
    try:
        mcube = shpmask.mask_cube(cube)
    except:
        print("Can't mask with shape! Masked over lon-lat box instead...")
        mcube = cube

    cubedata = []
    timedata = []

    # slicing cube over every ensemble run
    for subcube in mcube.slices_over('forecast_reference_time'):
        # extracting the time
        tcoord = subcube.coord('time')
        tdata = tcoord.units.num2date(tcoord.points.data)

        if 'max' in var:
            mcube = subcube.collapsed(['latitude', 'longitude'], iris.analysis.MAX)
        elif 'min' in var:
            mcube = subcube.collapsed(['latitude', 'longitude'], iris.analysis.MIN)
        else:
            raise NameError(f"No collapse method defined for {var}")

        cubedata.append(mcube.data)
        timedata.append(tdata)

    # Convert to Pandas Dataframe with unified time index
    s = list()
    [s.append(pd.Series(data=cubedata[i], index=timedata[i])) for i in range(np.shape(timedata)[0])]
    return pd.DataFrame(s).T


def get_era5_ts(indir, storm, shpmask):
    """
    Get ERA5 timeseries data of storm

    Arguments:
        indir (str): Path to ERA5 data
        storm (str): Name of storm

    Returns:
        Pandas dataframe
    """
    wfile = f"{indir}/era5.*.{storm}.global.hourly.nc"
    ucube = iris.load_cube(wfile, '10 metre U wind component')
    vcube = iris.load_cube(wfile, '10 metre V wind component')
    gube = iris.load_cube(wfile, '10 metre wind gust since previous post-processing')
    pcube = iris.load_cube(wfile, 'air_pressure_at_mean_sea_level')
    
    gube = gube.intersection(longitude=(75, 100), latitude=(10, 25))
    pcube = pcube.intersection(longitude=(75, 100), latitude=(10, 25))
    ucube = ucube.intersection(longitude=(75, 100), latitude=(10, 25))
    vcube = vcube.intersection(longitude=(75, 100), latitude=(10, 25))
    
    # Calculate wind speed
    ws_ifunc = iris.analysis.maths.IFunc(calc_wspd, ws_units_func)
    wcube = ws_ifunc(ucube, vcube, new_name='10 meter wind speed')

    try:
        mgcube = shpmask.mask_cube(gube)
        mpcube = shpmask.mask_cube(pcube)
        mwcube = shpmask.mask_cube(wcube)
    except:
        print("Can't mask with shape! Masked over lon-lat box instead...")
        mgcube = gube
        mpcube = pcube
        mwcube = wcube
    
    
    tcoord = mgcube.coord('time')
    tdata = tcoord.units.num2date(tcoord.points)

    gust = mgcube.collapsed(['latitude', 'longitude'], iris.analysis.MAX)
    pressure = mpcube.collapsed(['latitude', 'longitude'], iris.analysis.MIN)
    wind = mwcube.collapsed(['latitude', 'longitude'], iris.analysis.MAX)

    return pd.DataFrame(data={"GUST": mpsec_to_knots(gust.data), 
                              "PRES": pressure.data,
                              "WIND": mpsec_to_knots(wind.data)}, index=tdata)


def mpsec_to_knots(data):
    """
    function converts the wind units from meter per second to knots
    """
    return data * 1.945


def get_xtick_labels(stime, etime):
    """
    function returns the xtick labels
    """
    startday = stime[0].day
    endday = etime[-1, 8].astype(object).day + 3
    return range(startday, endday)



def get_ibtracs_uncertainity(wind):

    """
    Function calculates the uncertainity at each point of storm track, depending on the wind speed
    Source: https://www.ncdc.noaa.gov/ibtracs/pdf/IBTrACS_version4_Technical_Details.pdf
    Table 2 - Uncertainty of TC position based on TC intensity:
    wind < 60 kt                    ~ 30-40 km
    60 < wind < 100 kt              ~ 20-25 km
    wind > 100 kt                   ~ 10-15 km

    km to degree conversion:         1° = 111 km

    Arguments:
        Wind speed

    Returns:
        Uncertainity array
    """

    i = 0
    err = np.zeros(len(wind))

    for w in wind:

        if w < 60:
            err[i] = 0.36
        elif w >= 60 and w <= 100:
            err[i] = 0.22
        elif w > 100:
            err[i] = 0.1
        else: err[i] = 0

        i +=1

    return err


def calc_uncertainity_polygon(lon,lat,wspd):
    """
    Function creates the uncertainity polygon of storm track

    Arguments:
        Lon, Lat, Wind speed

    Returns:
        polygon
    """

    error = get_ibtracs_uncertainity(wspd)
    hours = np.arange(1,len(error)+1,1)


    # interpolate
    points_num = 100
    interp_hours = np.linspace(min(hours), max(hours), points_num)

    lon = interp1d(hours, lon, kind='linear')(interp_hours)
    lat = interp1d(hours, lat, kind='linear')(interp_hours)
    error = interp1d(hours, error, kind='linear')(interp_hours)



    # make polygon
    thetas = np.linspace(0, 2 * np.pi, 360)
    polygon_lon = lon[:,None] + error[:,None] * np.sin(thetas)
    polygon_lat = lat[:,None] + error[:,None] * np.cos(thetas)


    # circles
    ps = [Polygon(i) for i in np.dstack((polygon_lon, polygon_lat))]

    # list of convex hulls of subsequent circles
    n = range(len(ps)-1)
    convex_hulls = [MultiPolygon([ps[i], ps[i+1]]).convex_hull for i in n]

    # Final polygon
    polygons = cascaded_union(convex_hulls)

    lons, lats = polygons.exterior.xy

    return lons,lats




def geo_formatting(ax,extent,xticklst,yticklst):

    ax.set_xticks(xticklst, crs=ccrs.PlateCarree())
    ax.set_yticks(yticklst, crs=ccrs.PlateCarree())

    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()

    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    ax.coastlines(resolution='10m')
    ax.set_extent(extent)

    return None

def track_fitting(df):
    '''
    Function do polynimial fitting of order 3 and returns the fitted lat, lon
    '''
    df = df.reset_index()
    df = df.sort_values(by=[3], ascending=False)

    model = np.poly1d(np.polyfit(df[4], df[3], 3))
    lats = np.linspace(min(df[4]), max(df[4]), 30, endpoint = True)
    lons = model(lats)

    return lons,lats
