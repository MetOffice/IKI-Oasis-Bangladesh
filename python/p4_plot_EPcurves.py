'''
IKI Bangladesh: Plotting Exceedance probability curves 

Plotting exceedance probability (EP) curves based on fpgrand data.
Plots EP curves for 4 key locations + min and max EP curve.

Author: Hamish Steptoe (hamish.steptoe@metoffice.gov.uk)
Created: 
QA:
'''
import iris
import matplotlib.pyplot as plt
from user_vars import HCNC
import pandas as pd
import numpy as np

KEYLOCS = ["Cox's Bazar", "Chittagong", "Dhaka", "Comilla"]

# Load exceedance cube
cube = iris.load_cube(HCNC+'/fpgrand.exceedance_tmp.fg.T1Hmax.4p4km.nc')

gustthreshold = cube.coord('Gust speed classification threshold').points

# Load populated places
pop = pd.read_csv('sup/BGD_places.txt')

# Find Min and Max EP
lowestEP = cube.collapsed(['latitude', 'longitude'], iris.analysis.MIN).data
highestEP = cube.collapsed(['latitude', 'longitude'], iris.analysis.MAX).data

fig, ax = plt.subplots(figsize=(7, 4))

for index, row in pop.iterrows():
    # Get index of cube cell closest to pop locations
    ix = np.argmin(abs(cube.coord('longitude').points-row.LONGITUDE))
    iy = np.argmin(abs(cube.coord('latitude').points-row.LATITUDE))
    # Get exceedance probability
    ep = cube[iy, ix, :].data
    
    
    # Plotting
    if row.NAME in KEYLOCS:
        ax.semilogy(gustthreshold, ep, label=row.NAME, lw=2)
    else:
        ax.semilogy(gustthreshold, ep, color='black', alpha=0.3, lw=1)

ax.semilogy(gustthreshold, lowestEP, color='blue', ls='--', alpha=0.5)
ax.semilogy(gustthreshold, highestEP, color='blue', ls='--', alpha=0.5, label='Min/Max')
ax.grid(True, which='both', alpha=0.3)

plt.xlabel('Gust Speed (m/s)')
plt.ylabel('Storm Exceedance Probability')
plt.ylim([1e-3, 1])
plt.xlim([0, 80])
plt.legend(frameon=False)
plt.savefig('EPcurves.png', dpi=300)
plt.show()