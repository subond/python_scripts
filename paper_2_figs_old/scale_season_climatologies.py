"""
Load up the seasonal climatologies for year lengths >1, convert to scaled pentads (i.e. 10 days for double length, 20 for 4X, 40 for 8x)
"""

import xarray as xr
import sh
import numpy as np

#data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_2.000_unscaled.nc')

#data = data.rename({'xofyear':'pentad'})

#data.coords['xofyear'] = (data.pentad+1)//2

#print data.xofyear.values

#data = data.groupby('xofyear').mean(('pentad'))     

#filename = '/scratch/rg419/Data_moist/climatologies/sn_2.000.nc'
#data.to_netcdf(filename)



#data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_4.000_unscaled.nc')

#data = data.rename({'xofyear':'pentad'})

#data.coords['xofyear'] = (data.pentad+3)//4

#print data.xofyear.values

#data = data.groupby('xofyear').mean(('pentad'))     

#filename = '/scratch/rg419/Data_moist/climatologies/sn_4.000.nc'
#data.to_netcdf(filename)




data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_8.000_unscaled.nc')

data = data.rename({'xofyear':'pentad'})

data.coords['xofyear'] = (data.pentad+1)//2

print data.xofyear.values

data = data.groupby('xofyear').mean(('pentad'))     

filename = '/scratch/rg419/Data_moist/climatologies/sn_8.000.nc'
data.to_netcdf(filename)