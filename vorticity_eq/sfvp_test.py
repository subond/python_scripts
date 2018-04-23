"""
Test using windspharm for evaluating streamfunction and velocity potential for our data

"""
import cartopy.crs as ccrs
import matplotlib as mpl
mpl.rcParams['mathtext.default'] = 'regular'
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means

from windspharm.xarray import VectorWind
from windspharm.examples import example_data_path


# Read zonal and meridional wind components from file using the xarray module.
data = time_means('ap_1_partraw', [193,313], filename='atmos_pentad', timeav='month')

uwnd = data.ucomp
vwnd = data.vcomp

# Create a VectorWind instance to handle the computation of streamfunction and
# velocity potential.
w = VectorWind(uwnd, vwnd)

# Compute the streamfunction and velocity potential.
sf, vp = w.sfvp()
uchi, vchi, upsi, vpsi = w.helmholtz()

sf *= 1e-6
vp *= 1e-6

uchi[11,15,:,:].plot.contourf(x='lon',y='lat', extend='both')
plt.figure(2)
vchi[11,15,:,:].plot.contourf(x='lon',y='lat', extend='both')
plt.figure(3)
upsi[11,15,:,:].plot.contourf(x='lon',y='lat', extend='both')
plt.figure(4)
vpsi[11,15,:,:].plot.contourf(x='lon',y='lat', extend='both')
plt.show()

# Plot streamfunction.
#sf[11,15,:,:].plot.contourf(x='lon',y='lat',levels=np.arange(-120,121,20), extend='both')
# Plot velocity potential.
#plt.figure(2)
#vp[11,15,:,:].plot.contourf(x='lon',y='lat',levels=np.arange(-10,11,2), extend='both')
#plt.figure(3)
#uwnd[11,15,:,:].plot.contourf(x='lon',y='lat', extend='both')
# Plot velocity potential.
#plt.figure(4)
#vwnd[11,15,:,:].plot.contourf(x='lon',y='lat', extend='both')
#plt.show()

