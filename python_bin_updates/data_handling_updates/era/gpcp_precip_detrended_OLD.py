"""
Load in gpcp precip data. Detrend and create climatology

"""
import numpy as np
import xarray as xr
from scipy import signal

name_temp = '/scratch/rg419/obs_and_reanalysis/datafiles/gpcp_1dd_v1.2_p1d.%04d%02d.nc'
names = [name_temp % (m,n) for m in range( 1997, 2015) for n in range(1,13) ]

#read data into xarray 
data = xr.open_mfdataset( names, decode_times=False, chunks={'time': 30})

data.coords['day_since'] = (data.time * 0.) + np.arange(0., data.time.values.size)

data.swap_dims({'time':'day_since'}, inplace=True)

data.rename({'time':'day_of_yr'}, inplace=True)


# Locate leap years and discard from data
leap_list = np.argwhere(data.day_of_yr.values==365)
leap_list = [item for sublist in leap_list for item in sublist]

print 'start dropping'
data = data.drop(leap_list, dim='day_since')

# Add year number (year 0 is 1979)    
year_no =  np.repeat(np.arange(18.), 365)
year_no = (data.day_since * 0.) + year_no[0: data.day_since.values.size]
data.coords['year_no'] = year_no

# Make day_since a MultiIndex dimension, then unstack to get year_no and day_of_yr as separate array dimensions
data_rs = data['precip'].set_index(day_since = ['year_no', 'day_of_yr'])
data_rs = (data_rs.unstack('day_since'))
print 'Data reshaped'

# Set NaNs to zero and linearly detrend data along time axis
precip_detrend = signal.detrend(np.nan_to_num(data_rs), axis=data_rs.get_axis_num('year_no')) 
precip_detrend = xr.DataArray(precip_detrend, coords=data_rs.coords, dims=data_rs.dims)
# Add time mean back on
precip_detrend = precip_detrend + data_rs.mean('year_no')

# Convert day of year to pentad
pentad = precip_detrend.day_of_yr //5 + 1.

# Assign pentads as coordinate 'xofyear' to precip_detrend, and average to get climatology
precip_detrend.coords['xofyear'] = pentad
precip_clim = precip_detrend.groupby('xofyear').mean(('day_of_yr','year_no'))
precip_clim = precip_clim.transpose('xofyear', 'lat', 'lon')

ds = xr.Dataset({'precip_clim': (['xofyear', 'lat', 'lon'], precip_clim)},
                     coords={'xofyear': ('xofyear', precip_clim.xofyear),
                               'lat': ('lat', precip_clim.lat),
                               'lon': ('lon', precip_clim.lon)})
                              

# Write output climatology to netcdf
filename = '/scratch/rg419/obs_and_reanalysis/gpcp_detrendedpentads.nc'
ds.to_netcdf(filename)