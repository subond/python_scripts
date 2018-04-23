# 11/01/2018 Evaluate streamfunction, plot for JJA anomaly from time mean, and anom from zonal mean

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection



land = xr.open_dataset('/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc')

#(land.zsurf/1000).plot.contourf(x='lon', y='lat', levels = np.arange(0.,6.,0.5), add_colorbar=False, add_labels=False, extend='max', cmap='pink_r', alpha=0.8)
rcParams['figure.figsize'] = 12, 7.5
rcParams['font.size'] = 25
plot_dir = '/scratch/rg419/plots/era_wn2/'

(land.z[0,:,:]/9.8/1000).plot.contourf(x='longitude', y='latitude', levels = np.arange(0.,6.,0.5), add_colorbar=False, add_labels=False, extend='max', cmap='pink_r', alpha=0.8)
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=2)

#boxes = []
rect = Rectangle((70, 5), 30, 25, edgecolor='None', color='b', alpha=0.5)
ax = plt.gca()
ax.add_patch(rect)

rect = Rectangle((100, 20), 40, 20, edgecolor='None', color='b', alpha=0.5)
ax = plt.gca()
ax.add_patch(rect)

#rect = Rectangle((110, 30), 40, 15, edgecolor='None', color='k', alpha=0.5)
#ax = plt.gca()
#ax.add_patch(rect)

#rect = Rectangle((60, 30), 30, 15, edgecolor='None', color='k', alpha=0.5)
#ax = plt.gca()
#ax.add_patch(rect)

plt.xlim(50,170)
plt.ylim(-15,60)

plt.xlabel('Longitude')
plt.ylabel('Latitude')

plt.tight_layout()

plt.savefig(plot_dir + 'regions_map.pdf', format='pdf')
plt.close()