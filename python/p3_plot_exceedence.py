'''
IKI Bangladesh: Plotting Exceedance probability maps

Plotting spatial exceedance probabilities based on fpgrand data.

Author: HS
'''
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import iris
import iris.plot as iplt
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ascend import shape
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.colors import BoundaryNorm
from matplotlib.transforms import offset_copy
from mpl_toolkits.axes_grid1 import AxesGrid
from plotting import scale_bar
from user_vars import HCNC

# Load BGD outline
bgd = shape.load_shp('sup/BGD_shape.shp')
bgd_shape = bgd.unary_union()

# Load populated places
pop = pd.read_csv('sup/BGD_places.txt')

# Load BGD Level 1 admin regions
admin1 = shpreader.Reader('sup/BGD_admin1_lines_internal.shp')

# Load model data
gfp = iris.load_cube(HCNC + '/fpgrand.exceedance.fg.T1Hmax.4p4km.nc')

scs = iris.analysis.maths.exponentiate(gfp[:, :, 25], -1)
vscs = iris.analysis.maths.exponentiate(gfp[:, :, 33], -1)
sscs = iris.analysis.maths.exponentiate(gfp[:, :, 62], -1)

# Mask data to BGD shape
scs = bgd_shape.mask_cube(scs)
vscs.data.mask = scs.data.mask
sscs.data.mask = scs.data.mask

# Plotting
fig = plt.figure(figsize=(12, 8))
# Define grid of axes
axgr = AxesGrid(fig, 111, axes_class=(GeoAxes, dict(map_projection=ccrs.PlateCarree())),
                nrows_ncols=(1, 3), axes_pad=0.05,
                cbar_location='bottom', cbar_mode='single', cbar_pad=0.2, cbar_size='1%',
                label_mode='')  # note the empty label_mode

exceedenceprob = np.array([2, 5, 10, 20, 50, 100, 200])
levels = exceedenceprob
gustcmap = plt.get_cmap('plasma')
gustcmap.set_over('lightgrey', 1.0)
gustcmap.set_under('black', 1.0)
norm = BoundaryNorm(boundaries=levels, ncolors=gustcmap.N)

cf = iplt.pcolormesh(scs, axes=axgr[0], cmap=gustcmap, norm=norm)
iplt.pcolormesh(vscs, axes=axgr[1], cmap=gustcmap, norm=norm)
iplt.pcolormesh(sscs, axes=axgr[2], cmap=gustcmap, norm=norm)

# Plot place labels
for ax in axgr:
    ax.plot(pop['LONGITUDE'].values, pop['LATITUDE'].values, marker='o', fillstyle='none', markersize=5,
            linestyle='none', color='black',
            path_effects=[path_effects.Stroke(linewidth=2, foreground='white'), path_effects.Normal()])
    geodetic_transform = ccrs.PlateCarree()._as_mpl_transform(ax)
    for i, row, in pop.iterrows():
        if row['NAMEASCII'] in ['Dhaka', 'Saidpur', 'Jamalpur']:
            text = ax.text(row['LONGITUDE'], row['LATITUDE'], row['NAMEASCII'], ha='right', fontsize=8,
                           transform=offset_copy(geodetic_transform, units='dots', x=-4, y=4))
            text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='white', alpha=0.7),
                                   path_effects.Normal()])
        else:
            text = ax.text(row['LONGITUDE'], row['LATITUDE'], row['NAMEASCII'], ha='left', fontsize=8,
                           transform=offset_copy(geodetic_transform, units='dots', x=4, y=4))
            text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='white', alpha=0.7),
                                   path_effects.Normal()])
    # Add admin level 1 regions
    ax.add_geometries(admin1.geometries(), ccrs.PlateCarree(), fc='none', ec='white', lw=1, alpha=0.7)
    # Remove ax frames
    ax.background_patch.set_visible(False)
    ax.outline_patch.set_visible(False)
    # Add scale bar
    scale_bar(ax, (0.42, 0.1), 1_00)

axgr[0].text(0.5, 0.05, "Severe Cyclonic Storm\n(≥ 48 kt)\n(≥ 25 m/s)", ha='center', transform=axgr[0].transAxes,
             fontsize=10)
axgr[0].text(0.5, 0.05, "Very Severe Cyclonic Storm\n(≥ 64 kt)\n(≥ 33 m/s)", ha='center', transform=axgr[1].transAxes,
             fontsize=10)
axgr[1].text(0.5, 0.05, "Super Cyclonic Storm\n(≥ 120 kt)\n(≥ 62 m/s)", ha='center', transform=axgr[2].transAxes,
             fontsize=10)

cbar = fig.colorbar(cf, cax=axgr.cbar_axes[0], extend='both', orientation='horizontal')
cbar.ax.set_aspect('auto')
cbar.set_label("Event Exceedance Probability (%)")
cbar.ax.set_xticklabels(['50', '20', '10', '5', '2', '1', '0.5'])

plt.tight_layout()
plt.show()

plt.savefig('exceedance.png', dpi=300)
