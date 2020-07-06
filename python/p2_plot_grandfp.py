'''
IKI Bangladesh (MIOASI): Plot grand footprint

Plotting the 95th and 99th percentiles of the grand footprint data

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
bgd = shape.load_shp('BGD_shape.shp')
bgd_shape = bgd.unary_union()

# Load populated places
pop = pd.read_csv('BGD_places.txt')

# Load BGD Level 1 admin regions
admin1 = shpreader.Reader('BGD_admin1_lines_internal.shp')

# Load model data
gfp = iris.load_cube(HCNC + '/fpgrand.fg.T1Hmax.4p4km.nc', 'posterior credible intervals of wind_speed_of_gust')

q95 = gfp[:, :, -2]
q99 = gfp[:, :, -1]

# Mask data to BGD shape
q95 = bgd_shape.mask_cube(q95)
q99.data.mask = q95.data.mask

# Plotting
fig = plt.figure(figsize=(12, 8))
# Define grid of axes
axgr = AxesGrid(fig, 111, axes_class=(GeoAxes, dict(map_projection=ccrs.PlateCarree())),
                nrows_ncols=(1, 2), axes_pad=0.05,
                cbar_location='bottom', cbar_mode='single', cbar_pad=0.2, cbar_size='1%',
                label_mode='')  # note the empty label_mode

levels = np.arange(20, 65, 5)
gustcmap = plt.get_cmap('YlOrRd')
gustcmap.set_under('lightgrey', 1.0)
gustcmap.set_over('black', 1.0)
norm = BoundaryNorm(boundaries=levels, ncolors=gustcmap.N)

cf = iplt.pcolormesh(q95, axes=axgr[0], cmap=gustcmap, norm=norm)
iplt.pcolormesh(q99, axes=axgr[1], cmap=gustcmap, norm=norm)

# Plot Level 1 regions

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
    ax.add_geometries(admin1.geometries(), ccrs.PlateCarree(), fc='none', ec='white', lw=1)
    # Remove ax frames
    ax.background_patch.set_visible(False)
    ax.outline_patch.set_visible(False)
    # Add scale bar
    scale_bar(ax, (0.42, 0.17), 1_00)

axgr[0].text(0.5, 0.1, "95th percentile\n(1-in-20)", ha='center', transform=axgr[0].transAxes, fontsize=10)
axgr[1].text(0.5, 0.1, "99th percentile\n(1-in-100)", ha='center', transform=axgr[1].transAxes, fontsize=10)

cbar = fig.colorbar(cf, cax=axgr.cbar_axes[0], extend='both', orientation='horizontal')
cbar.set_label("Gust Speed (m/s)")


plt.tight_layout()
plt.show()

plt.savefig('grandfpq95q99.png', dpi=300)
