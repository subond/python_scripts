# Evaluate 'leading order' terms in vorticity equation: mean horizontal vel dotted with mean gradient of abs vort, and mean abs vort times mean vel divergence

from finite_difference import cfd
import matplotlib.pyplot as plt
import os
import xarray as xr
import numpy as np
from physics import gradients as gr


def vort_eq(run, month, filename='plev_daily', lev=150.):

    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
    name = name_temp % month 
    #read data into xarray 
    data = xr.open_dataset( name, decode_times=False)
        
    # Calculate planetary vorticity
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    #vor_rel = v_dx - u_dy
    
    dvordx = gr.ddx(vor)
    dvordy = gr.ddy(vor, vector=False)
    dvordp = gr.ddp(vor)
    
    # Calculate horizontal material derivative
    horiz_md = -1. * (data.ucomp * dvordx + data.vcomp * dvordy)
    
    # Calculate vertical material derivative
    vert_md = -1. * data.omega * dvordp
    
    # Now do the terms involving gradients of windspeed
    div = gr.ddx(data.ucomp) + gr.ddy(data.vcomp)
    stretching = -1. * vor * div
    tilting = gr.ddp(data.ucomp) * gr.ddy(data.omega, vector=False) - gr.ddp(data.vcomp) * gr.ddx(data.omega)
    
    total = horiz_md + vert_md + stretching + tilting
    
    ds = xr.Dataset({'horiz_md': (['time', 'pfull', 'lat', 'lon'], horiz_md),
                     'vert_md':  (['time', 'pfull', 'lat', 'lon'], vert_md),
                     'stretching':  (['time', 'pfull', 'lat', 'lon'], stretching),
                     'tilting':  (['time', 'pfull', 'lat', 'lon'], tilting),
                     'total':  (['time', 'pfull', 'lat', 'lon'], total)},
                     coords={'time': ('time', data.time),
                             'pfull': ('pfull', data.pfull),
                               'lat': ('lat', data.lat),
                               'lon': ('lon', data.lon)})
                              
    ds.coords['xofyear'] = np.mod( ds.time -1., 360.) //5 + 1.  
    dsout = ds.groupby('xofyear').mean(('time'))
    
    fileout = '/scratch/rg419/Data_moist/' + run + '/run%03d/vort_eq.nc'  
    fileout = fileout % month
    dsout.to_netcdf(path=fileout)
    
    print 'data written to ', fileout

    return dsout

for i in range(121,481):
    #vort_eq('ap_2', i)
    #vort_eq('full_qflux', i)
    vort_eq('flat_qflux', i)
    vort_eq('am_qflux', i)

#print 'Vorticity budget terms evaluated'
#data.horiz_md.plot()
#data.vert_md.plot()
#data.stretching.plot()
#data.tilting.plot()
#data.total.plot()
#plt.legend(['horiz', 'vert', 'stretching', 'tilting', 'residual'])
#plt.show()

