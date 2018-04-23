# Evaluate 'leading order' terms in vorticity equation: mean horizontal vel dotted with mean gradient of abs vort, and mean abs vort times mean vel divergence

from finite_difference import cfd
import matplotlib.pyplot as plt
import os
import xarray as xr
import numpy as np
from physics import gradients as gr


def vort_eq(run, months, filename='plev_daily', lev=150.):

    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    #read data into xarray 
    data = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    
    print 'files opened'
    
    # Calculate planetary vorticity
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    #vor_rel = v_dx - u_dy
    print 'absolute vorticity calculated'
    
    dvordx = gr.ddx(vor.sel(pfull=lev))
    dvordy = gr.ddy(vor.sel(pfull=lev), vector=False)
    dvordp = gr.ddp(vor)
    
    # Calculate horizontal material derivative
    horiz_md = -1. * (data.ucomp * dvordx + data.vcomp * dvordy)
    
    # Calculate vertical material derivative
    vert_md = -1. * data.omega.sel(pfull=lev) * dvordp.sel(pfull=lev)
    print 'advective terms done'
    
    # Now do the terms involving gradients of windspeed
    div = gr.ddx(data.ucomp.sel(pfull=lev)) + gr.ddy(data.vcomp.sel(pfull=lev))
    stretching = -1. * vor * div
    tilting = (gr.ddp(data.ucomp)).sel(pfull=lev) * gr.ddy(data.omega.sel(pfull=lev), vector=False) - (gr.ddp(data.vcomp)).sel(pfull=lev) * gr.ddx(data.omega.sel(pfull=lev))
    print 'stretching and tilting done'
    
    total = horiz_md + vert_md + stretching + tilting
    
    ds = xr.Dataset({'horiz_md': (['time', 'lat', 'lon'], horiz_md),
                     'vert_md':  (['time', 'lat', 'lon'], vert_md),
                     'stretching':  (['time', 'lat', 'lon'], stretching),
                     'tilting':  (['time', 'lat', 'lon'], tilting),
                     'total':  (['time', 'lat', 'lon'], total)},
                     coords={'time': ('time', data.time),
                             'pfull': ('pfull', data.pfull),
                               'lat': ('lat', data.lat),
                               'lon': ('lon', data.lon)})
                              
    ds.coords['xofyear'] = np.mod( ds.time -1., 360.) //5 + 1.  
    dsout = ds.groupby('xofyear').mean(('time'))
    
    GFDL_DATA = os.environ['GFDL_DATA']
    fileout = GFDL_DATA + 'climatologies/' + run + '_' + str(int(lev)) + '_vort_eq.nc'   
    dsout.to_netcdf(path=fileout)
    
    print 'data written to ', fileout

    return dsout

data = vort_eq('ap_2', [369,481])
data = vort_eq('full_qflux', [369,481])


#print 'Vorticity budget terms evaluated'
#data.horiz_md.plot()
#data.vert_md.plot()
#data.stretching.plot()
#data.tilting.plot()
#data.total.plot()
#plt.legend(['horiz', 'vert', 'stretching', 'tilting', 'residual'])
#plt.show()

