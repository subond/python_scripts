# Load in data for steady state simulations and save average of years 11-15 in climatology folder. 

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def load_steady_state(run, write_netcdf=True):
    
    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%04d/plev_pentad.nc'
    names = [name_temp % m for m in range( 61, 181)  ]
    #read data into xarray 
    data = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    
    data = data.mean('time')
    
    if write_netcdf:
        filename = '/scratch/rg419/Data_moist/climatologies/'+run+'.nc'
        data.to_netcdf(filename)
    

load_steady_state('sn_1.000_ss_250')
load_steady_state('sn_1.000_ss_260')
load_steady_state('sn_1.000_ss_270')

#load_steady_state('ss_94.000')    
#load_steady_state('ss_95.000')
#load_steady_state('ss_100.000')
#load_steady_state('ss_105.000')
#load_steady_state('ss_110.000')
#load_steady_state('ss_115.000')
