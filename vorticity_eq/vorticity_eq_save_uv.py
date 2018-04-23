# Save climatology of vorticity eq terms

from finite_difference import cfd
import matplotlib.pyplot as plt
import os
import xarray as xr
import numpy as np
from physics import gradients as gr


def save_vort_uv(run):
    #Load in dataset
    
    name = '/scratch/rg419/Data_moist/climatologies/' + run + '_vcomp_daily.nc'
    #read data into xarray 
    data_v = xr.open_dataset( name, decode_times=False)
    
    name = '/scratch/rg419/Data_moist/climatologies/' + run + '_ucomp150.0_daily.nc'
    #read data into xarray 
    data_u = xr.open_dataset( name, decode_times=False)
    
    print 'data loaded'
    
    data_v.coords['xofyear'] = np.mod( data_v.time -1., 360.) //5 + 1.  
    data_u.coords['xofyear'] = np.mod( data_u.time -1., 360.) //5 + 1.  
    
    dsout = data_v.groupby('xofyear').mean(('time'))
    dsout['ucomp'] = data_u.ucomp.groupby('xofyear').mean(('time'))
    
    #data_v.coords['pentad'] = data_v.coords['xofyear']
    #data_u.coords['pentad'] = data_u.coords['xofyear']
    
    #dsout = data_v.groupby('pentad').mean(('xofyear'))
    #dsout['ucomp'] = data_u.groupby('pentad').mean(('xofyear'))
    print dsout
    
    fileout = '/scratch/rg419/Data_moist/climatologies/vort_eq_uv' + run + '.nc'  
    dsout.to_netcdf(path=fileout, format = 'NETCDF3_CLASSIC')
    print 'data written to ', fileout


save_vort_uv('ap_2')
save_vort_uv('full_qflux')

