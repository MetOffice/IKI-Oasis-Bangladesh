'''
IKI Bangladesh (MIOASI): Storm tracks validation plotting

Plotting storm track density with IBTrACS tracks for IKI tropical cyclone events.

Note: Python 3 compatible only
Run using scitools/experimental-current (as of 14/5/20)

Author: ZM & HS
Created: 27/4/20
QA: HS
'''

import numpy as np
import matplotlib.pyplot as plt
import dataprocessing as dp
import cartopy.crs as ccrs
import glob
import pandas as pd
import shapely.geometry as sgeom
import rasterio
from rasterio import features
from ciid_tools import tc_tools
from ascend import shape
import cartopy.io.shapereader as shpreader
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
import matplotlib.patheffects as path_effects
from plotting import scale_bar
import pdb

IBTRACSFILE = 'sup/ibtracs.NI.list.v04r00.csv'
INDIR = 'storm_tracks/'

EVENTS = pd.read_csv('sup/events.csv', header=4)

# Load BGD outline
bgd = shape.load_shp('sup/NorthBoB_coast.shp')
bgd_shape = bgd.unary_union()

# Load populated places
pop = pd.read_csv('sup/BGD_places.txt')

# Load BGD Level 1 admin regions
admin1 = shpreader.Reader('sup/BGD_admin1_lines_internal.shp')
admin2 = shpreader.Reader('sup/WDBII_border_i_L1.shp')

# Load 4.4km geotiff to get rasterio dimensions
raster = rasterio.open('sup/4p4domain_30p0.tif')
# Get lon-lat coordinates - quirk of the function means I have to get them individually
llon, _ = raster.xy(np.zeros(raster.width), range(raster.width))
_, llat = raster.xy(range(raster.height), np.zeros(raster.height))
lons, lats = np.meshgrid(llon, llat)

# Set-up plotting
fig = plt.figure(figsize=(6.5, 7))
# Define grid of axes
axgr = AxesGrid(fig, 111, axes_class=(GeoAxes, dict(map_projection=ccrs.PlateCarree())),
                nrows_ncols=(4, 3), axes_pad=0.05, share_all=True,
                cbar_location='bottom', cbar_mode='single', cbar_pad=0.2, cbar_size='5%',
                label_mode='')  # note the empty label_mode
cmap = plt.get_cmap('viridis_r')
locs = np.arange(1, 10, 1)
norm = BoundaryNorm(boundaries=locs, ncolors=cmap.N)

# Sort in alphabetical order
EVENTS.sort_values(by='NAME', inplace=True)

for sn, STORM in enumerate(EVENTS.NAME.values):
    print(f'Plotting {STORM}...')

    # Rasterise tracks
    stracks = sorted(glob.glob(INDIR + f'/STITCHED_{STORM}_*_U-BJ924.DAT'))
    trackrast = np.zeros(raster.shape)
    for esmb in stracks:
        tetrack = tc_tools.TeFile(esmb, 'u-bj924')
        if not tetrack.tracks:
            break
        for track in tetrack.tracks:
            line = sgeom.LineString(list(zip(track.real_lon,
                                             track.real_lat)))
            # Convert Line to raster & add to raster array
            trackrast += features.rasterize([line], out_shape=raster.shape, transform=raster.transform,
                                           all_touched=True)
    
    # plot raster
    dx = raster.transform.a / 2.
    dy = raster.transform.e / 2.
    cf = axgr[sn].pcolormesh(lons - dx, lats - dy, np.ma.masked_less(trackrast, 1), snap=True, norm=norm, cmap=cmap,
                             zorder=5)
    
    # loading ibtracs data
    id = EVENTS['IBID'][EVENTS['NAME'] == STORM].values[0]
    ibtracs = dp.get_ibtracs_ts(IBTRACSFILE, id)
    ustracs = ibtracs[pd.notnull(ibtracs['USA_WIND'])]
    ndtracs = ibtracs[pd.notnull(ibtracs['NEWDELHI_WIND'])]
    

    # creating uncertainity polygons
    if not ustracs.empty:
        xs,ys = dp.calc_uncertainity_polygon(ustracs.USA_LON,ustracs.USA_LAT,ustracs.USA_WIND)
        axgr[sn].fill(xs, ys, fc='none', ec="#11dbdb", ls='--', lw=1, zorder=10)
        axgr[sn].plot(ustracs.USA_LON, ustracs.USA_LAT, color='white', lw=1, label='IBTrACS US',
                      path_effects=[path_effects.Stroke(linewidth=3, foreground='#11dbdb'), path_effects.Normal()],
                      zorder=10)
    if not ndtracs.empty:
        xd,yd = dp.calc_uncertainity_polygon(ndtracs.NEWDELHI_LON,ndtracs.NEWDELHI_LAT,ndtracs.NEWDELHI_WIND)
        axgr[sn].fill(xd, yd, fc='none', ec="#ff43ff", ls='--', lw=1, zorder=10)
        axgr[sn].plot(ndtracs.NEWDELHI_LON, ndtracs.NEWDELHI_LAT, color='white', lw=1, label='IBTrACS ND',
                      path_effects=[path_effects.Stroke(linewidth=3, foreground="#ff43ff"), path_effects.Normal()],
                      zorder=10)

    
    # Add coastline
    axgr[sn].add_geometries(bgd_shape.data, ccrs.PlateCarree(), fc='black', ec='none', lw=0, zorder=1, alpha=0.1)
    # Add admin level 1 regions
    axgr[sn].add_geometries(admin1.geometries(), ccrs.PlateCarree(), fc='none', ec='grey', lw=0.5, zorder=10)
    axgr[sn].add_geometries(admin2.geometries(), ccrs.PlateCarree(), fc='none', ec='grey', lw=0.5, zorder=10)
    
    # formatting the figure
    axgr[sn].text(0.1, 0.9, STORM, ha='left', transform=axgr[sn].transAxes, fontsize=8, weight='bold', zorder=20,
                  path_effects=[path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
    axgr[sn].set_extent([86.5, 94, 19.5, 25.0])
    # Remove ax frames
    axgr[sn].background_patch.set_visible(False)
    axgr[sn].outline_patch.set_visible(False)
    

# Add scale bar
scale_bar(axgr[-2], (0.2, 0.1), 1_00)
cbar = fig.colorbar(cf, cax=axgr.cbar_axes[0], orientation='horizontal', extend='max')
cbar.set_ticks(locs+.5)
cbar.set_ticklabels(locs)
cbar.ax.set_aspect('auto')
cbar.set_label("Storm Tracks per bin")
plt.tight_layout()
# plt.show()
plt.savefig('trackvalidation.png', dpi=600)
plt.close()
