"""Take model output and rewrite to netcdf"""

import numpy as np
from data_handling import copy_netcdf_attrs
import os
from netCDF4 import Dataset
import xarray as xr
from windspharm.xarray import VectorWind

GFDL_DATA = os.environ['GFDL_DATA']

filename = GFDL_DATA + 'full_qflux/run121/atmos_pentad.nc'
fileout = GFDL_DATA + 'full_qflux/run121/atmos_test.nc'

dsin= Dataset(filename, 'r', format='NETCDF3_CLASSIC')


data = xr.open_dataset(filename,decode_times=False)
uwnd = data.ucomp
vwnd = data.vcomp
# Create a VectorWind instance to handle the computation of streamfunction and
# velocity potential.
w = VectorWind(uwnd, vwnd)

# Compute the streamfunction and velocity potential.
streamfun, vel_pot = w.sfvp()


dsout= Dataset(fileout, 'w', format='NETCDF3_CLASSIC')

dsout = copy_netcdf_attrs(dsin, dsout)


sf_out = dsout.createVariable('streamfun', 'f4', ('time', 'pfull', 'lat', 'lon',))
sf_out.setncatts({k: dsout.variables['ps'].getncattr(k) for k in dsout.variables['ps'].ncattrs()})
sf_out.setncattr('long_name', 'streamfunction')
sf_out.setncattr('units', 'm**2 s**-1')
sf_out[:] = streamfun.load().data

vp_out = dsout.createVariable('vel_pot', 'f4', ('time', 'pfull', 'lat', 'lon',))
vp_out.setncatts({k: dsout.variables['ps'].getncattr(k) for k in dsout.variables['ps'].ncattrs()})
vp_out.setncattr('long_name', 'velocity potential')
vp_out.setncattr('units', 'm**2 s**-1')
vp_out[:] = vel_pot.load().data

dsout.close()

dsin.close()


from gfdl import Experiment

exp = Experiment('full_qflux')
#exp.runinterp(1,'atmos_pentad.nc','plev_pentad.nc', p_even=True)

exp.runinterp(121,'atmos_test.nc','plev_pentad_test.nc', p_even=True)
