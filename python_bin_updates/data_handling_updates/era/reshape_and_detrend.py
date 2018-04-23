"""
Load era data on single level
Find day of year
Reshape to year, day of year, lat, lon
Loop over day of year and pfull with detrend

"""
import numpy as np
import xarray as xr
from scipy import signal
import matplotlib.pyplot as plt
import time
import pandas as pd


def reshape_detrend(var, lev=None):
    
    # Open up data for a given variable
    if lev==None:
        filename = '/scratch/rg419/obs_and_reanalysis/era_' + var + '.nc'
    else:
        filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_' + var + '/era_' + var + '_' + str(lev) + '.nc'
    
    data = xr.open_dataset(filename, decode_times=False, chunks={'latitude': 100, 'longitude': 100})
    print 'data loaded'

  
    # Data starts in Jan 79, first get number of days since start
    day_since0179 = (data.time - np.min(data.time))//24 
    data.coords['daysince0179'] = day_since0179
    # Take average to get daily data from 6 hourly
    data_daily = data.groupby('daysince0179').mean('time')    
    print 'daily mean taken'
    
    # Use length of each year to assign a day of year to each day
    yrs = np.arange(1979,2018)
    yr_length = np.where( np.mod(yrs, 4.)==0., 366, 365)
    
    day_of_yr = np.arange(0.)
    for i in np.nditer(yr_length):
        day_of_yr = np.append(day_of_yr, np.arange(0.,i))
        
    # Convert to xarray coord, and shorten to only cover days in dataset 
    day_of_yr = (data_daily.daysince0179 * 0.) + day_of_yr[0: data_daily.daysince0179.values.size]
    data_daily.coords['day_of_yr'] = day_of_yr    
    
    # Locate leap years and discard from data
    leap_list = np.argwhere(day_of_yr.values==365)
    leap_list = [item for sublist in leap_list for item in sublist]
    print 'start dropping'
    data_daily = data_daily.drop(leap_list, dim='daysince0179')
    
    # Add year number (year 0 is 1979)    
    year_no =  np.repeat(np.arange(40.), 365)
    year_no = (data_daily.daysince0179 * 0.) + year_no[0: data_daily.daysince0179.values.size]
    data_daily.coords['year_no'] = year_no
    
    # Make daysince0179 a MultiIndex dimension, then unstack to get year_no and day_of_yr as separate array dimensions
    data_rs = data_daily[var].set_index(daysince0179 = ['year_no', 'day_of_yr'])
    
    # Select years 0 to 37, as 2017 is only a partial year
    if lev==None:
        data_rs = (data_rs.unstack('daysince0179'))[:,:,0:38,:]
    else:
        data_rs = (data_rs.unstack('daysince0179'))[:,:,:,0:38,:]
        print data_rs
    print 'Data reshaped'
            
    # Detrend data!
    data_detrend = signal.detrend(np.nan_to_num(data_rs), axis=data_rs.get_axis_num('year_no')) 
    data_detrend = xr.DataArray(data_detrend, coords=data_rs.coords, dims=data_rs.dims)
    # Add time mean back on
    data_detrend = data_detrend + data_rs.mean('year_no')
        
    #data_rs.sel(level=lev)[:,:,0,0].plot.contourf()
    #plt.figure(2)
    #data_detrend.sel(level=lev)[:,:,0,0].plot.contourf()
    #plt.figure(3)
    #data_rs.sel(level=lev)[:,:,30,0].plot.contourf()
    #plt.figure(4)
    #data_detrend.sel(level=lev)[:,:,30,0].plot.contourf()
    #plt.figure(5)
    #data_rs.sel(level=lev)[:,:,0,300].plot.contourf()
    #plt.figure(6)
    #data_detrend.sel(level=lev)[:,:,0,300].plot.contourf()
    #plt.show()
    print data_detrend    
    print 'Data detrended, saving'
    # Write detrended output to netcdf
    if lev==None:
        filename = '/scratch/rg419/obs_and_reanalysis/era_' + var + '_dt.nc'
    else:
        filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_' + var + '/era_' + var + '_' + str(lev) + '_dt.nc'
    
    if lev==None:
        ds = xr.Dataset({ var: (['lat', 'lon', 'year_no', 'day_of_yr'], data_detrend)},
                         coords={'lat': ('lat', data_detrend.latitude),
                                 'lon': ('lon', data_detrend.longitude),
                                 'year_no': ('year_no', data_detrend.year_no),
                                 'day_of_yr': ('day_of_yr', data_detrend.day_of_yr)})
    else:
        ds = xr.Dataset({ var + '_' + str(lev): (['level', 'lat', 'lon', 'year_no', 'day_of_yr'], data_detrend)},
                         coords={'level': ('level', data_detrend.level),
                                 'lat': ('lat', data_detrend.latitude),
                                 'lon': ('lon', data_detrend.longitude),
                                 'year_no': ('year_no', data_detrend.year_no),
                                 'day_of_yr': ('day_of_yr', data_detrend.day_of_yr)})
                          
    ds.to_netcdf(filename)


    
for var in ['u']:
    for lev in range(400,1001,50):
        t = time.time()
        reshape_detrend(var, lev=lev)
        elapsed = time.time() - t
        print var, lev, elapsed

