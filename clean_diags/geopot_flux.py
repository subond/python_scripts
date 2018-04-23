# Evaluate 'leading order' terms in vorticity equation: mean horizontal vel dotted with mean gradient of abs vort, and mean abs vort times mean vel divergence

from finite_difference import cfd
import matplotlib.pyplot as plt
import os
import xarray as xr
import numpy as np
from physics import gradients as gr


def geopot_flux(run, month, filename='plev_daily'):

    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
    name = name_temp % month 
    #read data into xarray 
    data = xr.open_dataset( name, decode_times=False)
    
    # Calculate geopotential flux
    vphi = data.vcomp * data.height * 9.8
    
    ds = xr.Dataset({'vphi': (['time', 'pfull', 'lat', 'lon'], vphi),
                     'vcomp':  (['time', 'pfull', 'lat', 'lon'], data.vcomp),
                     'phi':  (['time', 'pfull', 'lat', 'lon'], data.height*9.8)},
                     coords={'time': ('time', data.time),
                             'pfull': ('pfull', data.pfull),
                               'lat': ('lat', data.lat),
                               'lon': ('lon', data.lon)})
                              
    ds.coords['xofyear'] = np.mod( ds.time -1., 360.) //5 + 1.  
    dsout = ds.groupby('xofyear').mean(('time'))
    
    fileout = '/scratch/rg419/Data_moist/' + run + '/run%03d/vphi.nc'  
    fileout = fileout % month
    dsout.to_netcdf(path=fileout)
    
    print 'data written to ', fileout

    return dsout

for i in range(121,481):
    geopot_flux('full_qflux', i)
    


def save_geopot_flux(run, months):
    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/vphi.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    #read data into xarray 
    data = xr.open_mfdataset( names, decode_times=False)
    
    print 'data loaded'
    
    data.coords['pentad'] = data.coords['xofyear']
    
    dsout = data.groupby('pentad').mean(('xofyear'))
    
    fileout = '/scratch/rg419/Data_moist/climatologies/geopot_flux_' + run + '.nc'  
    dsout.to_netcdf(path=fileout, format = 'NETCDF3_CLASSIC')
    print 'data written to ', fileout

save_geopot_flux('full_qflux', [121,481])

