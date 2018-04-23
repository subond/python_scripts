"""
Load single level in era data. Take pentad means, detrend, and create climatology

"""
import numpy as np
import xarray as xr
from scipy import signal
import matplotlib.pyplot as plt


def load_era(var):
    
    # Open up data for a given variable
    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_'+var+'.nc', decode_times=False)
    
    print 'data loaded'
    
    # Data starts in Jan 79, first get number of days since start
    day_since0179 = (data.time - np.min(data.time))/24 
    data.coords['daysince0179'] = day_since0179
    # Take average to get daily data from 6 hourly
    data_daily = data.groupby('daysince0179').mean('time')
    
    print 'daily mean taken'
    
    # Set NaNs to zero and linearly detrend data along time axis
    data_detrend = signal.detrend(np.nan_to_num(data_daily[var]), axis=0) 
    data_detrend = xr.DataArray(data_detrend, coords=data_daily[var].coords, dims=data_daily[var].dims)
    # Add time mean back on
    data_detrend = data_detrend + data_daily[var].mean('daysince0179')
    
    print 'detrending complete'
    
    # Data covers years 1979-2017, work out how long these years are
    yrs = np.arange(1979,2018)
    yr_length = np.where( np.mod(yrs, 4.)==0., 366, 365)
    
    # Use length of each year to assign a day of year to each day
    day_of_yr = np.arange(0.)
    for i in np.nditer(yr_length):
        day_of_yr = np.append(day_of_yr, np.arange(0.,i))
    
    # Convert to xarray coord, and shortern to only cover days in dataset (2017 has only Jan and Feb)
    day_of_yr = (data_detrend.daysince0179 * 0.) + day_of_yr[0: data_detrend.daysince0179.values.size]
    # Add as coord in data_detrend 
    data_detrend.coords['day_of_yr'] = day_of_yr
    
    # Convert day of year to pentad
    pentad = data_detrend.day_of_yr //5 + 1.
    
    # Find locations of any pentad 74s (leap years)
    leap_list = np.argwhere(pentad.values==74.)
    # Create an array of [-1,0,0,0,0,-1,....] to add from Feb-Dec of leap years so that extra day is in pentad 12
    pentad_correction = np.zeros([306,])
    pentad_correction[0::5] = -1
    # Loop over leap years identified and add correction
    for leap in leap_list:
        pentad[leap[0]-305:leap[0]+1] = pentad[leap[0]-305:leap[0]+1] + pentad_correction
    
    # Assign a pentad number to each day. Pentad 12 has 6 days in a leap year
    data_detrend.coords['xofyear'] = pentad
    # Take pentad mean to get climatology
    data_clim = data_detrend.groupby('xofyear').mean('daysince0179')
    
    print 'climatology calculated, saving...'
    
    ds = xr.Dataset({var: (['xofyear', 'lat', 'lon'], data_clim)},
                         coords={'xofyear': ('xofyear', data_clim.xofyear),
                                   'lat': ('lat', data_clim.latitude),
                                   'lon': ('lon', data_clim.longitude)})
                              

    # Write output climatology to netcdf
    filename = '/scratch/rg419/obs_and_reanalysis/era_clim_' + var + '.nc'
    ds.to_netcdf(filename)
    

load_era('z')
load_era('d')