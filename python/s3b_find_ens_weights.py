'''
IKI Bangladesh (MIOASI): S3b Find weights for ens footprints

Find weights for use in s3_sample_storm_ensemble.R based on assessing
the similarity of the fpens ensemble members to the ERA5 equivalent footprint.

Similarity is based on the DISO metric, which combines r (correlation coefficient),
AE (absolute error, measuring any persistent bias) and RMSE (root mean square error,
averaged magnitude of the deviation) to summarise model performance. Each metric
is given equal weight such that the DISO is the Euclidean distance between the obs
and the 3 indices in 3D space defined by r, AE and RMSE.

Author: HS
Created: 10/12/19
'''

import glob
import iris
import re
from iris.analysis.stats import pearsonr
from iris.analysis.maths import abs as cabs
import csv
from user_vars import HCNC


def diso(obs, model):
    '''DISO (Distance between Indices of Simulation and Observation)
    Index defined by Hu et al. (2018, doi:10.1002/joc.5972) designed to to describe
    the overall performances of different models against the observed field quantitatively.

    It merges r (correlation coefficient), AE (absolute error, measuring any persistent bias)
    and RMSE (root mean square error, averaged magnitude of the deviation) to summarise model
    performance. Each metric is given equal weight such that the DISO is the Euclidean distance
    between the obs and the 3 indices in 3D space defined by r, AE and RMSE.

    [If this doesn't make sense, go read the paper - it has a nice figure]

    In DISO space, diso(obs) = 0.  diso(model) = 0 indicates that the model exactly
    matches the obs, according to metrics measured by RMSE, AE and r.

    Args:
        obs (iris.cube.Cube): Cube of base values (to which you are comparing)
        model (iris.cube.Cube): Cube of model values

    Returns:
        (rr, nae, nrmse) (tupple): Stats used to calculate the DISO
        diso (iris.cube.Cube): Cube of DISO values
    '''
    # Check input
    if not isinstance(obs, iris.cube.Cube) or not isinstance(model, iris.cube.Cube):
        raise TypeError('Args should be Iris cubes')
    # Check bounds
    if not obs.coord(axis='x').has_bounds():
        obs.coord(axis='x').guess_bounds()
        obs.coord(axis='y').guess_bounds()
    if not model.coord(axis='x').has_bounds():
        model.coord(axis='x').guess_bounds()
        model.coord(axis='y').guess_bounds()
    # Correlation coefficient
    rr = pearsonr(obs, model)
    rr.convert_units('1')
    # Difference
    diff = model - obs
    # Absolute obs mean
    gridweights = iris.analysis.cartography.area_weights(obs)
    # TODO: Generalise collapsing of cube dimensions, handle rotated pole and time?
    aobs_bar = cabs(obs.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=gridweights))
    # Normalised Absolute Error
    nae = diff.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=gridweights) / aobs_bar
    nae.rename('Normalised Absolute Error')
    # Normalised RMSE
    nrmse = rmse(model, obs) / aobs_bar
    nrmse.rename('Normalised Root Mean Square Error')
    # Distance between Indices of Simulation and Observation
    diso = ((rr - 1) ** 2 + nae ** 2 + nrmse ** 2) ** 0.5
    diso.rename('Distance between Indices of Simulation and Observation')
    return (rr, nae, nrmse), diso

def rmse(obs, model):
    '''Iris Root Mean Square

    Iris Root Mean Sqaure (RMS) function.  When applied to the model - obs difference,
    this gives you the Root Mean Square Error (RMSE).

    Args:
        obs (iris.cube.Cube): Cube of base values (to which you are comparing)
        model (iris.cube.Cube): Cube of model values

    Returns:
        Iris cube

    '''
    # Check input
    if not isinstance(obs, iris.cube.Cube) or not isinstance(model, iris.cube.Cube):
        raise TypeError('Args should be Iris cubes')
    diff = diff_cube(model, obs)
    if not diff.coord(axis='x').has_bounds():
        diff.coord(axis='x').guess_bounds()
        diff.coord(axis='y').guess_bounds()
    w = iris.analysis.cartography.area_weights(diff)
    rms = diff.collapsed(['latitude', 'longitude'], iris.analysis.RMS, weights=w)
    rms.rename('Root Mean Square')
    return rms

def diff_cube(model, obs, obs_name=None, ensmem=None):
    '''Create a difference cube = model - obs,
    rename cube appropriately.
    '''
    diff = model - obs
    diff.standard_name = model.standard_name
    if obs_name is not None and ensmem is not None:
        diff.long_name = '{}_{}{}_{}_diff'.format(model.standard_name, user_vars.MODEL_NAME,
                                                  ensmem, obs_name)
    return diff

# Set boundary for data intersection
LAT = (20.503502, 27.483002)
LON = (87.5555, 92.942)
# Set ERA5 var name
VAR = "10 metre wind gust since previous post-processing"

# Get fpens file names
fpensfiles = glob.glob(HCNC + '/fpens.fg.T1Hmax.*.4p4km.nc')
# Open CSV file
with open('fpens_r.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerow(['Name', 'Date', 'r'])
    # Loop through fpens files
    for file in fpensfiles:
        name = re.search(r'.*T1Hmax.(.*?).4p4.*', file).group(1)
        # Load downscaled ensemble
        ens = iris.load_cube(file).intersection(longitude=LON, latitude=LAT)
        # Load ERA5
        era5 = iris.load_cube(f'fp.era5.*.{name}.nc', VAR).intersection(longitude=LON, latitude=LAT)
        # Add coordinate system
        era5.coord(axis='x').coord_system = ens.coord(axis='x').coord_system
        era5.coord(axis='y').coord_system = ens.coord(axis='y').coord_system
        # Regrid ERA5 onto ens grid
        era5r = iris.util.squeeze(era5.regrid(ens[0], iris.analysis.Nearest()))
        for e in ens.slices_over('forecast_reference_time'):
            frt = e.coord('forecast_reference_time')
            frts = frt.units.num2date(frt.points)[0].strftime('%Y-%m-%d %H')
            # metric = 1 - diso(e, era5r)[1].data
            metric = pearsonr(e, era5r).data
            print(f'{name} {frts}: {metric:.4f}')
            csvwriter.writerow([name, frts, metric])
    






