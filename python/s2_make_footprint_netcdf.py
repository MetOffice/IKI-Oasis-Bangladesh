"""
IKI Project Plotting: Make footprint ensemble file

For each downscaled storm, make a netCDF file of the the footprint of each
ensemble member

Author: HS
Created: 16/7/19
"""
import argparse
import glob
import iris
from user_vars import EVENTS, HCNC, RES, VARS, DOMAINS

# Parse other processing options from the command line
parser = argparse.ArgumentParser(description="Specify file options")
parser.add_argument(
    "-v", "--variable", choices=list(VARS.keys()), default="fg.T1Hmax", help="Specify variable to process"
)
parser.add_argument("-r", "--resolution", choices=RES, default="4p4", help="Specify pp file resolution")
parser.add_argument(
    "-e",
    "--event",
    choices=list(EVENTS.keys()),
    nargs="*",
    default=list(EVENTS.keys()),
    help="Specify event (storm name)",
)
pargs = parser.parse_args()

# Define the statistic used for making the footprints
FPSTAT = iris.analysis.MAX

for e in pargs.event:
    # Load data
    file = glob.glob(f"{HCNC}/{pargs.variable}*{e.upper()}.{pargs.resolution}km.nc")
    cube = iris.load_cube(file)

    # Make common domain size for 1.5km simulations
    if pargs.resolution == "1p5":
        cube = cube.intersection(
            longitude=(DOMAINS["BGD"][0], DOMAINS["BGD"][1]), latitude=(DOMAINS["BGD"][2], DOMAINS["BGD"][3])
        )
    # To collapse across the forecast period, we need to remove the time dimensions
    try:
        cube.remove_coord("time")
    except iris.exceptions.CoordinateNotFoundError:
        pass
    # Make footprint over forecast_period, ie. over each ensemble member
    fp = cube.collapsed("forecast_period", VARS[pargs.variable])
    # Save footprint cub
    outfile = f"fpens.{pargs.variable}.{e.upper()}.{pargs.resolution}km.nc"
    print(f"Saving {outfile} to {HCNC}...")
    iris.save(fp, f"{HCNC}/{outfile}")
