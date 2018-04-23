# Evaluate 'leading order' terms in vorticity equation: mean horizontal vel dotted with mean gradient of abs vort, and mean abs vort times mean vel divergence

import matplotlib.pyplot as plt
import os
import xarray as xr
import numpy as np
from pylab import rcParams
import sh

rcParams['figure.figsize'] = 10, 6
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/steady_state_runs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
    
def sst_plot(qlat, qamp):
    
    #read data into xarray 
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/qflux_' + qlat + '_' + qamp + '.nc')
    
    t_anom = data.t_surf - data.t_surf.mean('lon')
    p_anom = data.ps - data.ps.mean('lon')
    u_anom = data.ucomp.sel(pfull=850.) - data.ucomp.mean('lon').sel(pfull=850.)
    v_anom = data.vcomp.sel(pfull=850.) - data.vcomp.mean('lon').sel(pfull=850.)
    
    t_anom.plot.contourf(x='lon', y='lat', levels=np.arange(-8.,8.1,1.), add_labels=False)
    #data_q = xr.open_dataset('/scratch/rg419/python_scripts/qflux_anoms/qflux_150_lat_' + qlat + '_a_' + qamp + '.nc')
    #data_q.ocean_qflux.plot(x='lon', y='lat', levels=np.arange(0.,300.,100.), colors='k')
    plt.quiver(data.ucomp.lon[::3], data.ucomp.lat[::1], u_anom[::1,::3], v_anom[::1,::3])#, scale=500.,headwidth=5)
    plt.contour(p_anom.lon, p_anom.lat, p_anom, add_labels=False, colors='0.7', levels=np.arange(-300.,301.,50.), linewidths=2)
    plt.xlim([50., 250.])
    plt.ylim([-45.,45.])
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    figname = 'sst_' + qlat + '_' + qamp + '.pdf'
    plt.savefig(plot_dir + figname)
    plt.close()


    
def w_plot(qlat, qamp):
    
    #read data into xarray 
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/qflux_' + qlat + '_' + qamp + '.nc')
    
    t_anom = data.t_surf - data.t_surf.mean('lon')
    w_anom = (data.omega - data.omega.mean('lon')).sel(pfull=500.)*86400./100. # Convert to hPa/day
    #w_anom.plot.contourf(x='lon', y='lat')
    #plt.show()
    
    u_anom = data.ucomp.sel(pfull=850.) - data.ucomp.mean('lon').sel(pfull=850.)
    v_anom = data.vcomp.sel(pfull=850.) - data.vcomp.mean('lon').sel(pfull=850.)
    
    t_anom.plot.contourf(x='lon', y='lat', levels=np.arange(-8.,8.1,1.), add_labels=False)
    #data_q = xr.open_dataset('/scratch/rg419/python_scripts/qflux_anoms/qflux_150_lat_' + qlat + '_a_' + qamp + '.nc')
    #data_q.ocean_qflux.plot(x='lon', y='lat', levels=np.arange(0.,300.,100.), colors='k')
    plt.quiver(data.ucomp.lon[::3], data.ucomp.lat[::1], u_anom[::1,::3], v_anom[::1,::3])#, scale=500.,headwidth=5)
    plt.contour(w_anom.lon, w_anom.lat, w_anom, add_labels=False, colors='0.7', levels=np.arange(-180.,181.,20.), linewidths=2)
    plt.xlim([50., 250.])
    plt.ylim([-45.,45.])
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    figname = 'w_' + qlat + '_' + qamp + '.pdf'
    plt.savefig(plot_dir + figname)
    plt.close()


w_plot('0', '300')
w_plot('5', '300')
w_plot('10', '300')
w_plot('15', '300')
w_plot('20', '300')
w_plot('25', '300')




