'''
IKI Bangladesh: Plotting decision map

Plot map of optional decision based on loss function.

Author: HS
'''
import iris
import numpy as np
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import iris.plot as iplt
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from ascend import shape
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.colors import BoundaryNorm, from_levels_and_colors
from matplotlib.transforms import offset_copy
from mpl_toolkits.axes_grid1 import AxesGrid
from plotting import scale_bar
import pandas as pd


def getprobs(posterior):
    """
    Find probabilities given defined thresholds.
    In this case, 4 hard-coded threshold levels.
    Args:
        posterior: posterior predictive distribution (iris.Cube)

    Returns: np.array of threshold probabilities

    """
    
    prob1 = np.mean(posterior < 13.9, axis=2)
    prob2 = np.mean((posterior >= 13.9) & (posterior < 16.9), axis=2)
    prob3 = np.mean((posterior >= 16.9) & (posterior < 24.8), axis=2)
    prob4 = np.mean((posterior >= 24.8), axis=2)
    # NB. The sum of these probabilities = 1
    return np.array([prob1, prob2, prob3, prob4])


# Define relative losses per threshold per action
LOSSFUN = np.array([[0, 5, 15, 20],  # grey column
                    [50, 10, 20, 25],  # yellow column
                    [80, 50, 25, 30],  # orange column
                    [100, 100, 80, 40]])  # red column

# Based on Economou et al. (2016) paper
# LOSSFUN = np.array([[0, 4, 12, 25],  # grey column
#                     [70, 38, 28, 25],  # yellow column
#                     [88, 46, 31, 25],  # orange column
#                     [100, 52, 34, 25]])  # red column

predscube = iris.load_cube('r/rdata/posterior_predictions.nc')
# Get probability of occurences for each of the intervals defined in getprobs()
# NB.
probs = getprobs(predscube.data)

# Find loss associated with each action
# ie. matrix multiplication -> in 2D this is LOSSFUN (4x4) * probs (1x4) -> loss per action (1x4)
# in python: np.matmul(LOSSFUN.T, probs) OR (LOSSFUN.T*probs).sum(axis=1) OR np.einsum('ij,i->j',LOSSFUN, probs)
# to broadcast in 3D for four probs levels ie. (4, lon, lat) needs some matrix voodoo help from np.einsum()
# see https://stackoverflow.com/questions/26089893/understanding-numpys-einsum
# NB. Very important to get the dimensions the right way round
actionloss = np.einsum('ij,ikl->jkl', LOSSFUN, probs)
# then find action with lowest loss (ie. along axis=0) with np.argmin(),
# ie. which action minimises the loss
warn = np.argmin(actionloss, axis=0)

# Load BGD outline
bgd = shape.load_shp('sup/BGD_shape.shp')
bgd_shape = bgd.unary_union()

# Mask data
# Set-up dummy cube for warn data
warncube = predscube[:,:,0]
warncube.data = warn
warncube.rename('warnings')
warncube = bgd_shape.mask_cube(warncube)

# Load populated places
pop = pd.read_csv('sup/BGD_places.txt')

# Load BGD Level 1 admin regions
admin1 = shpreader.Reader('sup/BGD_admin1_lines_internal.shp')


# Plotting
fig = plt.figure(figsize=(6, 8))
# Define grid of axes
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

levels = np.array([0, 1, 2, 3, 4])
cmap, norm = from_levels_and_colors(levels, ['grey', 'wheat', 'orange', 'red'])

cf = iplt.pcolormesh(warncube, axes=ax, cmap=cmap, norm=norm)

# Plot place labels

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
# scale_bar(ax, (0.42, 0.13), 1_00)

cbar = fig.colorbar(cf, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04,
                    ticks=levels+0.5)
cbar.outline.set_visible(False)
cbar.solids.set_edgecolor('white')
cbar.solids.set_linewidth(5)
cbar.ax.set_aspect('auto')
cbar.ax.set_xlabel('Warning Level')
cbar.ax.set_xticklabels(['None', 'Alert/Warning', 'Disaster', 'Great Danger'])

fig.tight_layout()
# plt.show()

plt.savefig('warnings.png', dpi=300)