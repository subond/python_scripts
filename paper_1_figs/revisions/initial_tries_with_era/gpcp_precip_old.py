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
#data = xr.open_mfdataset( names, chunks={'time': 30})

# Set NaNs to zero and linearly detrend data along time axis
precip_detrend = signal.detrend(np.nan_to_num(data.precip), axis=0) 
precip_detrend = xr.DataArray(precip_detrend, coords=data.precip.coords, dims=data.precip.dims)
# Add time mean back on
precip_detrend = precip_detrend + data.precip.mean('time')

# Convert day of year to pentad
pentad = data.time //5 + 1.
# Find locations of any pentad 74s (leap years)
leap_list = np.argwhere(pentad.values==74.)
# Create an array of [-1,0,0,0,0,-1,....] to add from Feb-Dec of leap years so that extra day is in pentad 12
pentad_correction = np.zeros([306,])
pentad_correction[0::5] = -1
# Loop over leap years identified and add correction
for leap in leap_list:
    pentad[leap[0]-305:leap[0]+1] = pentad[leap[0]-305:leap[0]+1] + pentad_correction


# Assign pentads as coordinate 'xofyear' to precip_detrend, and average to get climatology
precip_detrend.coords['xofyear'] = pentad
precip_clim = precip_detrend.groupby('xofyear').mean('time')

ds = xr.Dataset({'precip_clim': (['xofyear', 'lat', 'lon'], precip_clim)},
                     coords={'xofyear': ('xofyear', precip_clim.xofyear),
                               'lat': ('lat', precip_clim.lat),
                               'lon': ('lon', precip_clim.lon)})
                              

# Write output climatology to netcdf
filename = '/scratch/rg419/obs_and_reanalysis/gpcp_detrendedpentads.nc'
ds.to_netcdf(filename)