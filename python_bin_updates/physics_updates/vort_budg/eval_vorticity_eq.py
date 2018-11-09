# Evaluate 'leading order' terms in vorticity equation: mean horizontal vel dotted with mean gradient of abs vort, and mean abs vort times mean vel divergence

from finite_difference import cfd
import matplotlib.pyplot as plt
import os
import xarray as xr
import numpy as np
from data_handling_updates import gradients as gr


def vort_eq(run, month, filename='plev_daily', period_fac=1., rot_fac=1., do_ss=False, day_fac=1.):
    '''Evaluate terms in the vorticity budget using daily data from simulations. RG 2/11/17
       Inputs: run = run name
               month = month to evaluate
               filename = name of file to look for, default is pressure interpolated daily data
               period_fac = ratio of the orbital period to Earth's
               rot_fac = ratio of the rotation rate to Earth's
               do_ss = evaluate for a steady state simulation, do not average daily data, average all saved data together
       '''
    
    #Load in dataset
    try:
        name_temp = '/disca/share/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
        name = name_temp % month 
        #read data into xarr
        data = xr.open_dataset( name, decode_times=False)
    except:
        name_temp = '/disca/share/rg419/Data_moist/' + run + '/run%04d/'+filename+'.nc'
        name = name_temp % month 
        #read data into xarr
        data = xr.open_dataset( name, decode_times=False)
        
    # Calculate planetary vorticity
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    
    # Take gradients of vorticity
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
                     'ucomp': (['time', 'pfull', 'lat', 'lon'], data.ucomp),
                     'vcomp': (['time', 'pfull', 'lat', 'lon'], data.vcomp),
                     'total':  (['time', 'pfull', 'lat', 'lon'], total)},
                     coords={'time': ('time', data.time),
                             'pfull': ('pfull', data.pfull),
                               'lat': ('lat', data.lat),
                               'lon': ('lon', data.lon)})
    
    # If a steady state run, leave time dimension in to average out later. Otherwise take pentad means for climatology
    if do_ss:
        dsout = ds
    else:
        ds.coords['xofyear'] = np.mod( ds.time/day_fac -1., 360.*period_fac) //5 + 1.  
        dsout = ds.groupby('xofyear').mean(('time'))
    
    #try:
    #    fileout = '/disca/share/rg419/Data_moist/' + run + '/run%03d/vort_eq.nc'  
    #    fileout = fileout % month
    #    dsout.to_netcdf(path=fileout)
    #except:
    fileout = '/disca/share/rg419/Data_moist/' + run + '/run%04d/vort_eq.nc'  
    fileout = fileout % month
    dsout.to_netcdf(path=fileout)
    
    print ('data written to ', fileout)

    return dsout
    
                              

def save_vort_eq(run, months, do_ss=False):
    '''Load in vorticity files saved by the above, and produce a climatology. RG 2/11/2017
       Inputs: run: name of run
               months: start and end month to load
               do_ss: average data from a steady state simulation'''
    
    #Load in dataset
    try:
        name_temp = '/disca/share/rg419/Data_moist/' + run + '/run%03d/vort_eq.nc'
        names = [name_temp % m for m in range( months[0], months[1])  ]
        #read data into xarray 
        data = xr.open_mfdataset( names, decode_times=False)
    except:
        name_temp = '/disca/share/rg419/Data_moist/' + run + '/run%04d/vort_eq.nc'
        names = [name_temp % m for m in range( months[0], months[1])  ]
        #read data into xarray 
        data = xr.open_mfdataset( names, decode_times=False)
    
    print ('data loaded')
    
    # If a steady state run, average over all times. Otherwise add a second pentad coordinate and use it to create a clmatology 
    if do_ss:
        dsout = data.mean('time')
    else:
        data.coords['pentad'] = data.coords['xofyear']
        dsout = data.groupby('pentad').mean(('xofyear'))
    
    fileout = '/disca/share/rg419/Data_moist/climatologies/vort_eq_' + run + '.nc'  
    dsout.to_netcdf(path=fileout, format = 'NETCDF3_CLASSIC')
    print ('data written to ', fileout)


if __name__ == "__main__":
    
    #rot_fac=[0.75,1.,1.25,1.5]
    for run in ['mld_20', 'mld_15']:   
        k=0
        for i in range(121,481):
            vort_eq(run, i)
        save_vort_eq(run, [121,481])
        k=k+1
        