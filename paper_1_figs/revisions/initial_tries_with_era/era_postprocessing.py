"""
Load era data
Find day of year
Reshape to year, day of year, pfull, lat, lon
Loop over day of year and pfull with detrend
Take climatology

"""
import numpy as np
import xarray as xr
from scipy import signal
import matplotlib.pyplot as plt
import time
import pandas as pd



def day_av(var, twodim = True):
    
    # Open up data for a given variable
    t = time.time()
    
    if twodim:
        data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_'+var+'.nc', decode_times=False, chunks={'latitude': 100, 'longitude': 100})
    else:
        names = ['/scratch/rg419/obs_and_reanalysis/era_' + var + '_1979_1988.nc',
                 '/scratch/rg419/obs_and_reanalysis/era_' + var + '_1989_1998.nc',
                 '/scratch/rg419/obs_and_reanalysis/era_' + var + '_1999_2008.nc',
                 '/scratch/rg419/obs_and_reanalysis/era_' + var + '_2009_2017.nc']
        data = xr.open_mfdataset(names, decode_times=False, chunks={'latitude': 10, 'longitude': 10})
        #data = data.sel(level=range(50, 1001, 50))
    
    print 'data loaded'
  
    # Data starts in Jan 79, first get number of days since start
    day_since0179 = (data.time - np.min(data.time))/24 
    data.coords['daysince0179'] = day_since0179
    # Take average to get daily data from 6 hourly
    data_daily = data.groupby('daysince0179').mean('time')
    
    print 'daily mean taken'

    yrs = np.arange(1979,2018)
    yr_length = np.where( np.mod(yrs, 4.)==0., 366, 365)
    
    # Use length of each year to assign a day of year to each day
    day_of_yr = np.arange(0.)
    for i in np.nditer(yr_length):
        day_of_yr = np.append(day_of_yr, np.arange(0.,i))
        
    # Convert to xarray coord, and shortern to only cover days in dataset (2017 has only Jan and Feb)
    day_of_yr = (data_daily.daysince0179 * 0.) + day_of_yr[0: data_daily.daysince0179.values.size]
    # Add as coord in data_detrend 
    data_daily.coords['day_of_yr'] = day_of_yr    
    
    leap_list = np.argwhere(day_of_yr.values==365)
    leap_list = [item for sublist in leap_list for item in sublist]
    
    print 'start dropping'
    data_daily = data_daily.drop(leap_list, dim='daysince0179')
    
    year_no =  np.repeat(np.arange(40.), 365)
    year_no = (data_daily.daysince0179 * 0.) + year_no[0: data_daily.daysince0179.values.size]
    data_daily.coords['year_no'] = year_no
    
    #print 'set index'
    test = data_daily[var].set_index(daysince0179 = ['year_no', 'day_of_yr'])
    print test
    
    test2 = test.unstack('daysince0179')
    
    print test2
    
    elapsed = time.time() - t
    print elapsed
    
    for lev in np.nditer(data_daily.level.values):
        print lev
        ds = xr.Dataset({ var + '_' + str(lev): (['lat', 'lon', 'year_no', 'day_of_yr'], test2.sel(level=lev))},
                             coords={'lat': ('lat', test2.latitude),
                                     'lon': ('lon', test2.longitude),
                                     'year_no': ('year_no', test2.year_no),
                                     'day_of_yr': ('day_of_yr', test2.day_of_yr)})
                              

        # Write output climatology to netcdf
        filename = '/scratch/rg419/obs_and_reanalysis/'+ var + '_' + str(lev) + '.nc'
        ds.to_netcdf(filename)
        
    
    #data_detrend = signal.detrend(np.nan_to_num(test2), axis=3) 
    #data_detrend = xr.DataArray(data_detrend, coords=test2.coords, dims=test2.dims)
    # Add time mean back on
    #data_detrend = data_detrend + test2.mean('day_of_yr')
    
    #print data_detrend


data = day_av('vo', twodim=False)
