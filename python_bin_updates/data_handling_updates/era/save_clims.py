"""
Load detrended era data and save climatologies of required parts
"""
import numpy as np
import xarray as xr
from scipy import signal
import matplotlib.pyplot as plt
import time
import pandas as pd


def one_var_clim(var, lev=None):
    # Save climatology of just one variable
    
    # Open up data for a given variable
    if lev==None:
        filename = '/scratch/rg419/obs_and_reanalysis/era_' + var + '_dt.nc'
    else:
        filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_' + var + '/era_' + var + '_' + str(lev) + '_dt.nc'
    
    data = xr.open_dataset(filename, decode_times=False, chunks={'lat': 100, 'lon': 100})
    print 'data loaded'
    
    data = data.mean('year_no').sel(level=lev)
  
    # Write climatology to netcdf
    if lev==None:
        filename = '/scratch/rg419/obs_and_reanalysis/era_' + var + '_clim.nc'
        ds = xr.Dataset({ var: (['lat', 'lon', 'day_of_yr'], data[var])},
                         coords={'lat': ('lat', data.lat),
                                 'lon': ('lon', data.lon),
                                 'day_of_yr': ('day_of_yr', data.day_of_yr)})
    else:
        filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_' + var + '/era_' + var + '_' + str(lev) + '_clim.nc'
        ds = xr.Dataset({ var + '_' + str(lev): (['lat', 'lon', 'day_of_yr'], data[var + '_' + str(lev)])},
                         coords={'lat': ('lat', data.lat),
                                 'lon': ('lon', data.lon),
                                 'day_of_yr': ('day_of_yr', data.day_of_yr)})
                                 
                          
    ds.to_netcdf(filename)



def prod_clim(var1, var2, lev1=None, lev2=None):
    # Save climatology of product of variables
    
    # Open up data for a given variable
    if lev1==None:
        filename1 = '/scratch/rg419/obs_and_reanalysis/era_' + var1 + '_dt.nc'
    else:
        filename1 = '/scratch/rg419/obs_and_reanalysis/sep_levs_' + var1 + '/era_' + var1 + '_' + str(lev1) + '_dt.nc'
        
    if lev2==None:
        filename2 = '/scratch/rg419/obs_and_reanalysis/era_' + var2 + '_dt.nc'
    else:
        filename2 = '/scratch/rg419/obs_and_reanalysis/sep_levs_' + var2 + '/era_' + var2 + '_' + str(lev2) + '_dt.nc'
    
    print filename1
    print filename2
    
    data1 = xr.open_dataset(filename1, decode_times=False, chunks={'lat': 100, 'lon': 100})
    data2 = xr.open_dataset(filename2, decode_times=False, chunks={'lat': 100, 'lon': 100})
    print 'data loaded'
  
    # Write climatology to netcdf
    if (lev1==None and lev2==None):
        data = (data1[var1] * data2[var2]).mean('year_no')
        
        filename = '/scratch/rg419/obs_and_reanalysis/era_' + var1 + var2 + '_clim.nc'
        ds = xr.Dataset({ var1 + var2 : (['lat', 'lon', 'day_of_yr'], data)},
                         coords={'lat': ('lat', data.lat),
                                 'lon': ('lon', data.lon),
                                 'day_of_yr': ('day_of_yr', data.day_of_yr)})
    elif (lev1==None and lev2!=None):
        data = (data1[var1] * data2[var2 + '_' + str(lev2)]).mean('year_no')
        
        filename = '/scratch/rg419/obs_and_reanalysis/era_' + var1 + var2 + '_' + str(lev2) + '_clim.nc'
        ds = xr.Dataset({ var1 + var2 + '_' + str(lev2): (['lat', 'lon', 'day_of_yr'], data)},
                         coords={'lat': ('lat', data.lat),
                                 'lon': ('lon', data.lon),
                                 'day_of_yr': ('day_of_yr', data.day_of_yr)})
    
    elif (lev1!=None and lev2==None):
        data = (data1[var1 + '_' + str(lev1)] * data2[var2]).mean('year_no')
        
        filename = '/scratch/rg419/obs_and_reanalysis/era_' + var1 + var2 + '_' + str(lev1) + '_clim.nc'
        ds = xr.Dataset({ var1 + var2 + '_' + str(lev1): (['lat', 'lon', 'day_of_yr'], data)},
                         coords={'lat': ('lat', data.lat),
                                 'lon': ('lon', data.lon),
                                 'day_of_yr': ('day_of_yr', data.day_of_yr)})
    
    else:
        data = (data1[var1 + '_' + str(lev1)] * data2[var2 + '_' + str(lev2)]).mean('year_no')
        
        filename = '/scratch/rg419/obs_and_reanalysis/era_' + var1 + var2 + '_' + str(lev1) + '_clim.nc'
        ds = xr.Dataset({ var1 + var2 + '_' + str(lev1): (['lat', 'lon', 'day_of_yr'], data)},
                         coords={'lat': ('lat', data.lat),
                                 'lon': ('lon', data.lon),
                                 'day_of_yr': ('day_of_yr', data.day_of_yr)})
                                 
                          
    ds.to_netcdf(filename)
    
    

#prod_clim('vo', 'd', lev1=150)
#prod_clim('u', 'u', lev1=150, lev2=150)
#prod_clim('u', 'v', lev1=150)
#for lev in range(50,350,50):
#    prod_clim('u', 'w', lev1=lev, lev2=lev)
    

#for var in ['z','d','v']:
#    t = time.time()
#    one_var_clim(var)
#    elapsed = time.time() - t
#    print var, elapsed
    
for var in ['u']:
    for lev in range(350,1001,50):
        t = time.time()
        one_var_clim(var, lev=lev)
        elapsed = time.time() - t
        print var, lev, elapsed
