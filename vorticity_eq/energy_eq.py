# Evaluate 'leading order' terms in vorticity equation: mean horizontal vel dotted with mean gradient of abs vort, and mean abs vort times mean vel divergence

from finite_difference import cfd
import matplotlib.pyplot as plt
import os
import xarray as xr
import numpy as np
from physics import gradients as gr


def energy_eq(run, month, filename='plev_daily', lev=150.):

    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
    name = name_temp % month 
    #read data into xarray 
    data = xr.open_dataset( name, decode_times=False)
    
    # Calculate planetary vorticity
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    uchi = data.ucomp - data.upsi
    vchi = data.vcomp - data.vpsi
    
    pe_to_ke = -1. * (uchi * gr.ddx(data.height*9.8) + vchi * gr.ddy(data.height*9.8, vector=False))
        
    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    
    conv_1 = vor * (data.upsi*vchi - uchi*data.vpsi)
    
    conv_2 = -1. * omega * (data.upsi * gr.ddp(uchi) + data.vpsi * gr.ddp(vchi))
    
    conv_3 = -1. * (gr.ddx(uchi) +gr.ddy(vchi)) * (data.upsi**2. + data.vpsi**2.)/2.
    
    ds = xr.Dataset({'pe_to_ke': (['time', 'pfull', 'lat', 'lon'], pe_to_ke),
                     'conv_1':  (['time', 'pfull', 'lat', 'lon'], conv_1),
                     'conv_2':  (['time', 'pfull', 'lat', 'lon'], conv_2),
                     'conv_3':  (['time', 'pfull', 'lat', 'lon'], conv_3),
                     'upsi':  (['time', 'pfull', 'lat', 'lon'], data.upsi),
                     'vpsi':  (['time', 'pfull', 'lat', 'lon'], data.vpsi),
                     'uchi':  (['time', 'pfull', 'lat', 'lon'], uchi),
                     'vchi':  (['time', 'pfull', 'lat', 'lon'], vchi)},
                     coords={'time': ('time', data.time),
                             'pfull': ('pfull', data.pfull),
                               'lat': ('lat', data.lat),
                               'lon': ('lon', data.lon)})
                              
    ds.coords['xofyear'] = np.mod( ds.time -1., 360.) //5 + 1.  
    dsout = ds.groupby('xofyear').mean(('time'))
    
    fileout = '/scratch/rg419/Data_moist/' + run + '/run%03d/energy_eq.nc'  
    fileout = fileout % month
    dsout.to_netcdf(path=fileout)
    
    print 'data written to ', fileout

    return dsout

for i in range(121,481):
    energy_eq('full_qflux', i)
    


def save_energy_eq(run, months):
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

save_energy_eq('full_qflux', [121,481])

