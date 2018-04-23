"""
Load in horizontal streamfunction and plot

"""
import os
import numpy as np
import xarray as xr
import sh
from pylab import rcParams
from data_handling import add_vars_netcdf


GFDL_DATA = os.environ['GFDL_DATA']
filename_temp = GFDL_DATA + 'full_qflux/run%03d/atmos_pentad.nc'
filenames = [filename_temp % m for m in range(121, 481)]

for filename in filenames:
    print filename
    add_vars_netcdf(filename, vp=False, ups=False, vps=False)



