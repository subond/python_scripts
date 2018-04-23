# Have a stab at evaluating the divergent and non divergent winds speeds from div and omega

from physics import mass_streamfunction
from data_handling import time_means
from finite_difference import cfd
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np


#Assume, whether rightly or wrongly, that these are 0 at the poles... Then divergent v can be calculated using dv/dy = D - du/dx, and nondiv u using du/dy = dv/dx - vort. Can then multiply by grid spacing and integrate to get winds, and subtract from total to get other part of wind... 

def wind_components(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1.):
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    dudx = xr.DataArray( cfd(data.ucomp.values,data.lon*np.pi/180.,3),   [('xofyear', data.xofyear ), ('pfull', data.pfull ), ('lat', data.lat), ('lon', data.lon )])
    dudx = dudx #/coslat/a Don't bother with geometric factor as will just be multiplied by in integration
    
    dvdx = xr.DataArray( cfd(data.vcomp.values,data.lon*np.pi/180.,3),   [('xofyear', data.xofyear ), ('pfull', data.pfull ), ('lat', data.lat), ('lon', data.lon )])
    dvdx = dvdx #/coslat/a
    
    dvdy_div = data.div - dudx
    dudy_vort = dvdx - data.vor
    
    dlat=xr.DataArray(data.latb.diff('latb').values*np.pi/180, coords=[('lat', data.lat)])
    
    v_div  = np.cumsum(dvdy_div  * dlat, axis=2)/coslat
    u_vort = np.cumsum(dudy_vort * dlat, axis=2)/coslat
    
    u_div = data.ucomp - u_vort
    v_vort = data.vcomp - v_div
    
    plt.figure(1)
    v_vort[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-10,11,1))
    plt.figure(2)
    v_div[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-10,11,1))
    plt.figure(3)
    data.vcomp[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-10,11,1))
    plt.figure(4)
    u_vort[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-80,81,10))
    plt.figure(5)
    u_div[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-80,81,10))    
    plt.figure(6)
    data.ucomp[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-80,81,10))
    plt.show()
    
wind_components('ap_1_partraw', [193,313])





def wind_comp_uv(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1.):
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    dudy = xr.DataArray( cfd((data.ucomp*coslat).values,data.lon*np.pi/180.,3),   [('xofyear', data.xofyear ), ('pfull', data.pfull ), ('lat', data.lat), ('lon', data.lon )])
    dvdy = xr.DataArray( cfd((data.vcomp*coslat).values,data.lon*np.pi/180.,3),   [('xofyear', data.xofyear ), ('pfull', data.pfull ), ('lat', data.lat), ('lon', data.lon )])
    
    dlat=xr.DataArray(data.latb.diff('latb').values*np.pi/180, coords=[('lat', data.lat)])
    
    v_div  = np.cumsum(dvdy * dlat, axis=2)/coslat
    u_vort = np.cumsum(dudy * dlat, axis=2)/coslat
    
    u_div = data.ucomp - u_vort
    v_vort = data.vcomp - v_div
    
    plt.figure(1)
    v_vort[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-10,11,1))
    plt.figure(2)
    v_div[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-10,11,1))
    plt.figure(3)
    data.vcomp[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-10,11,1))
    plt.figure(4)
    u_vort[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-80,81,10))
    plt.figure(5)
    u_div[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-80,81,10))    
    plt.figure(6)
    data.ucomp[40,15,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-80,81,10))
    plt.show()
    
wind_comp_uv('ap_1_partraw', [193,313])

