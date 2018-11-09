"""
Convert climatologies of individual levels to one netcdf
"""
import numpy as np
import xarray as xr
from scipy import signal
import matplotlib.pyplot as plt
import time
import pandas as pd


def var_clim(var):
    
    var_array = np.zeros([365,20,73,144])
    i=0
    for lev in range(50, 1050, 50):
        name = '/scratch/rg419/obs_and_reanalysis/sep_levs_' + var + '/era_' + var + '_' + str(lev) + '_clim.nc' 
        data = xr.open_dataset(name, decode_times=False)
        var_array[:,i,:,:] = np.moveaxis(data[var+'_'+str(lev)].values, [0,1,2], [1,2,0])
        i=i+1
    
    ds = xr.Dataset({ var : (['day_of_yr', 'pfull', 'lat', 'lon'], var_array)},
                     coords={'day_of_yr': ('day_of_yr', data.day_of_yr),
                             'pfull': ('pfull', np.arange(50.,1050.,50.)),
                             'lat': ('lat', data.lat),
                             'lon': ('lon', data.lon)})

    ds.to_netcdf('/scratch/rg419/obs_and_reanalysis/era_' + var + '_clim_alllevs.nc' )


    
for var in ['u','v']:
    var_clim(var)