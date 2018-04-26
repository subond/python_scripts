"""
Save a climatology for sn8 run, needs 3 climatologies making and averaging first as there's so much data 26/02/2018

"""
import xarray as xr
import numpy as np
from data_handling_updates import time_means

#name_out = '/scratch/rg419/Data_moist/climatologies/sn_8.000_1.nc'
#test1=time_means('sn_8.000', [193,1153], filename='plev_pentad', timeav='pentad', write_netcdf=True, period_fac=8., name_out=name_out)
#name_out = '/scratch/rg419/Data_moist/climatologies/sn_8.000_2.nc'
#test2=time_means('sn_8.000', [1153,2113], filename='plev_pentad', timeav='pentad', write_netcdf=True, period_fac=8., name_out=name_out)
#name_out = '/scratch/rg419/Data_moist/climatologies/sn_8.000_3.nc'
#test3=time_means('sn_8.000', [2113,3073], filename='plev_pentad', timeav='pentad', write_netcdf=True, period_fac=8., name_out=name_out)

test1 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_8.000_1.nc')
test2 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_8.000_2.nc')
test3 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_8.000_3.nc')


data = (test1 + test2 + test3)/3.
name_out = '/scratch/rg419/Data_moist/climatologies/sn_8.000.nc'
data.to_netcdf(name_out)