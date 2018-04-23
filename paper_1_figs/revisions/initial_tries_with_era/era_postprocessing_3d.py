"""
Load era data. Take pentad means, detrend, and create climatology

"""
import numpy as np
import xarray as xr
from scipy import signal
import matplotlib.pyplot as plt
import time


def day_av(var, twodim = True):
    
    # Open up data for a given variable
    t = time.time()

    if twodim:
        data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_'+var+'.nc', decode_times=False, chunks={'time': 360})
    else:
        names = ['/scratch/rg419/obs_and_reanalysis/era_' + var + '_1979_1988.nc',
                 '/scratch/rg419/obs_and_reanalysis/era_' + var + '_1989_1998.nc',
                 '/scratch/rg419/obs_and_reanalysis/era_' + var + '_1999_2008.nc',
                 '/scratch/rg419/obs_and_reanalysis/era_' + var + '_2009_2017.nc']
        data = xr.open_mfdataset(names, decode_times=False, chunks={'latitude': 10, 'longitude': 10})
    
    print 'data loaded'
    elapsed = time.time() - t
    print elapsed
    
    t = time.time()
    
    # Data starts in Jan 79, first get number of days since start
    day_since0179 = (data.time - np.min(data.time))/24 
    data.coords['daysince0179'] = day_since0179
    # Take average to get daily data from 6 hourly
    data_daily = data.groupby('daysince0179').mean('time')
    
    print 'daily mean taken'
    elapsed = time.time() - t
    print elapsed
    
    return data_daily
    #filename = '/scratch/rg419/obs_and_reanalysis/era_daily_' + var + '.nc'
    #data_daily.to_netcdf(filename)
        

def load_and_detrend(data, twodim = True):
    
    # Open up data for a given variable
    #data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_daily_'+var+'.nc', decode_times=False, chunks={'latitude': 100, 'longitude': 100})
    
    #print 'data loaded'
    
    t = time.time()
    # Set NaNs to zero and linearly detrend data along time axis
    data_detrend = signal.detrend(np.nan_to_num(data[var]), axis=0) 
    data_detrend = xr.DataArray(data_detrend, coords=data[var].coords, dims=data[var].dims)
    # Add time mean back on
    data_detrend = data_detrend + data[var].mean('daysince0179')
    
    print 'detrending complete'
    elapsed = time.time() - t
    print elapsed
    
    #filename = '/scratch/rg419/obs_and_reanalysis/era_detrend_' + var + '.nc'
    #data_detrend.to_netcdf(filename)
    
    return data_detrend


def calculate_clim(data_detrend, twodim = True):
    #Take detrended time series and calculate a climatology
    
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
    
    if twodim:
        ds = xr.Dataset({var: (['xofyear', 'lat', 'lon'], data_clim)},
                         coords={'xofyear': ('xofyear', data_clim.xofyear),
                                   'lat': ('lat', data_clim.latitude),
                                   'lon': ('lon', data_clim.longitude)})
    else:
        ds = xr.Dataset({var: (['xofyear', 'pfull', 'lat', 'lon'], data_clim)},
                         coords={'xofyear': ('xofyear', data_clim.xofyear),
                                 'pfull': ('pfull', data_clim.level),
                                   'lat': ('lat', data_clim.latitude),
                                   'lon': ('lon', data_clim.longitude)})
                              

    # Write output climatology to netcdf
    filename = '/scratch/rg419/obs_and_reanalysis/era_clim_' + var + '.nc'
    ds.to_netcdf(filename)
    
data = day_av('vo', twodim=False)
vo_detrend = load_and_detrend(data, twodim=False)
calculate_clim(vo_detrend, twodim=False)

data = day_av('w', twodim=False)
w_detrend = load_and_detrend(data, twodim=False)
calculate_clim(w_detrend, twodim=False)
