# Look at progression of peak of omega throughout the year

from physics import mass_streamfunction
from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

a=6376.0e3
omega = 7.2921150e-5

data_ap303 = time_means('sn_3.000', [145,469], filename='atmos_pentad', timeav='pentad',period_fac=3.)
data_ap101 = time_means('aquaplanet_10m', [109,121], filename='atmos_pentad', timeav='pentad',period_fac=1.)

data_ap101.omega[:,27,:,:].mean('lon').plot.contourf(x='xofyear', y='lat',add_labels=False)
plt.ylim(-45,45)
plt.xlabel('Pentad')
plt.ylabel('Latitude')
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/w_101.png')
plt.close()

data_ap303.omega[:,27,:,:].mean('lon').plot.contourf(x='xofyear', y='lat',add_labels=False)
plt.ylim(-45,45)
plt.xlabel('Pentad')
plt.ylabel('Latitude')
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/w_303.png')
plt.close()
