# Evaluate 'leading order' terms in vorticity equation: mean horizontal vel dotted with mean gradient of abs vort, and mean abs vort times mean vel divergence

from finite_difference import cfd
import matplotlib.pyplot as plt
import os
import xarray as xr
import numpy as np
from physics import gradients as gr


def save_vort_eq(run, months):
    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/energy_eq.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    #read data into xarray 
    data = xr.open_mfdataset( names, decode_times=False)
    
    print 'data loaded'
    
    data.coords['pentad'] = data.coords['xofyear']
    
    dsout = data.groupby('pentad').mean(('xofyear'))
    
    fileout = '/scratch/rg419/Data_moist/climatologies/energy_eq_' + run + '.nc'  
    dsout.to_netcdf(path=fileout, format = 'NETCDF3_CLASSIC')
    print 'data written to ', fileout


save_vort_eq('ap_2', [121,481])
#save_vort_eq('full_qflux', [121,481])

