'''
IKI Project Plotting: Domain and downscaling nest plot

Author: HS
Created: 21/6/19
'''
import matplotlib.pyplot as plt
import numpy as np
import cartopy
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from ascend import shape
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def xmin(key):
    return CENTLON-(DELTA[key]*(NUMPTS[key][1]/2)) + OFFSET[key][1]

def xmax(key):
    return CENTLON+(DELTA[key]*(NUMPTS[key][1]/2)) + OFFSET[key][1]

def ymin(key):
    return CENTLAT-(DELTA[key]*(NUMPTS[key][0]/2)) + OFFSET[key][0]

def ymax(key):
    return CENTLAT+(DELTA[key]*(NUMPTS[key][0]/2)) + OFFSET[key][0]

CENTLAT, CENTLON = [24, 89]
RIM = 13
# Region Names
RGNAME = {'1': ['4.5km', 'red'],
          '2': ['1.4km', 'blue'],
          '3': ['ERA5', 'black']}
# Region Grid size (degrees)
DELTA = {'1': 0.0405,
         '2': 0.0135,
         '3': 0.1}
# Number of Grid Points (y, x)
NUMPTS = {'1': [816, 836],
          '2': [544, 426],
          '3': [360, 360]}
# Offset in degrees (y, x)
OFFSET = {'1': [0,0],
          '2': [0, 1.2555],
          '3': [0,0]}

# Load Validation area shapefile
val = shape.load_shp('/project/ciid/projects/IKI/ArcGIS/ValidationArea2.shp')
val_shape = val.unary_union()

# Construct boundaries
bounds = {}
for key in RGNAME.keys():
    bounds[key] = {'withrim': np.array([[xmin(key), xmin(key), xmax(key), xmax(key), xmin(key)],
                               [ymin(key), ymax(key), ymax(key), ymin(key), ymin(key)]]),
                   'norim': np.array([[xmin(key)+DELTA[key]*RIM, xmin(key)+DELTA[key]*RIM, xmax(key)-DELTA[key]*RIM,
                                       xmax(key)-DELTA[key]*RIM, xmin(key)+DELTA[key]*RIM],
                             [ymin(key)+DELTA[key]*RIM, ymax(key)-DELTA[key]*RIM, ymax(key)-DELTA[key]*RIM,
                              ymin(key)+DELTA[key]*RIM, ymin(key)+DELTA[key]*RIM]])
                   }

ax = plt.subplot(projection=ccrs.PlateCarree())
# ax.stock_img()
land_50m = cartopy.feature.NaturalEarthFeature('physical', 'land', '50m',
                                               edgecolor='face',
                                               facecolor=cartopy.feature.COLORS['land'])
# pdb.set_trace()
ax.add_feature(land_50m)
ax.coastlines(resolution='50m', lw=0.5, alpha=0.5)
ax.add_feature(cartopy.feature.BORDERS, lw=0.5, alpha=0.5)

for key in ['1','2']:
    ax.plot(bounds[key]['norim'][0], bounds[key]['norim'][1], lw=2, color=RGNAME[key][1])
    # ax.fill(x, y, color='coral', alpha=0.2)
# Plot ERA5
ax.plot(bounds['3']['withrim'][0], bounds['3']['withrim'][1], lw=2, color=RGNAME['3'][1])

# Plot validation domain
val_shape.plot(ax=ax, facecolor='None', edgecolor='green', lw=2)

# Plot ticks
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', alpha=0.8, linestyle=(0, (5, 10)))
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 8, 'color': 'gray'}
gl.ylabel_style = {'size': 8, 'color': 'gray'}

# Plot LSM
# plt.pcolormesh(lsm.coord('longitude').contiguous_bounds(),
#                lsm.coord('latitude').contiguous_bounds(),
#                lsm.data, transform=cartopy.crs.PlateCarree(),
#                vmax=1, cmap='OrRd_r', zorder=100)

# ax.set_extent([-30, 45, 30, 76.5])
plt.tight_layout(pad=3)
plt.show()
plt.savefig('domain.eps', dpi=300)