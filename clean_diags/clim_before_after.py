"""
Plot up 150 and 850 hPa u and v, and precip and surface temperature before and after onset

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means
from physics import mass_streamfunction
import sh
from finite_difference import cfd
from pylab import rcParams


    
def clim_ba(run, before_after, land=False):
    
    rcParams['figure.figsize'] = 15, 10
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if land:
        land_data = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')
        
    
    t_surf_before = data.t_surf[before_after[0]:before_after[0]+4,:,:].mean('xofyear')    
    t_surf_after = data.t_surf[before_after[1]:before_after[1]+4,:,:].mean('xofyear')
    precip_before = 86400.*(data.condensation_rain + data.convection_rain)[before_after[0]:before_after[0]+4,:,:].mean('xofyear')
    precip_after = 86400.*(data.condensation_rain + data.convection_rain)[before_after[1]:before_after[1]+4,:,:].mean('xofyear')
    u_before = data.ucomp[before_after[0]:before_after[0]+4,:,:,:].mean('xofyear')
    u_after = data.ucomp[before_after[1]:before_after[1]+4,:,:,:].mean('xofyear')
    v_before = data.vcomp[before_after[0]:before_after[0]+4,:,:,:].mean('xofyear')
    v_after = data.vcomp[before_after[1]:before_after[1]+4,:,:,:].mean('xofyear')
    

    # Four subplots
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    #First plot
    f1 = precip_before.plot.contourf(x='lon', y='lat', ax=ax1, extend = 'both', levels = np.arange(2.,14.1,2.), add_colorbar=False, add_labels=False)
    if land:
        land_data.land_mask.plot.contour(x='lon', y='lat', ax=ax1, colors='w', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False)
    ax1.quiver(data.lon[21:54:3], data.lat[10:55:1], u_before[17,10:55:1,21:54:3], v_before[17,10:55:1,21:54:3], headlength=5, headwidth=5)    
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)
    ax1.grid(True,linestyle=':')
    
    #Second plot
    f2 = precip_after.plot.contourf(x='lon', y='lat', ax=ax2, extend = 'both', levels = np.arange(2.,14.1,2.), add_colorbar=False, add_labels=False)
    if land:
        land_data.land_mask.plot.contour(x='lon', y='lat', ax=ax2, colors='w', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False)
    ax2.quiver(data.lon[21:54:3], data.lat[10:55:1], u_after[17,10:55:1,21:54:3], v_after[17,10:55:1,21:54:3], headlength=5, headwidth=5)
    ax2.grid(True,linestyle=':')
    
    #Third plot
    f3 = t_surf_before.plot.contourf(x='lon', y='lat', ax=ax3, extend = 'both', levels=range(250,321,10), add_colorbar=False, add_labels=False)
    if land:
        land_data.land_mask.plot.contour(x='lon', y='lat', ax=ax3, colors='w', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False)
    ax3.quiver(data.lon[21:54:3], data.lat[10:55:1], u_before[3,10:55:1,21:54:3], v_before[3,10:55:1,21:54:3], headlength=5, headwidth=5)
    ax3.grid(True,linestyle=':')
    ax3.set_ylabel('Latitude')
    ax3.set_xlabel('Longitude')
    ax3.set_ylim(-60,60)
    ax3.set_xlim(60,150)
    ax3.set_xticks(range(60,151,15))
    
    #Fourth plot
    f4 = t_surf_after.plot.contourf(x='lon', y='lat', ax=ax4, extend = 'both', levels=range(250,321,10), add_colorbar=False, add_labels=False)
    if land:
        land_data.land_mask.plot.contour(x='lon', y='lat', ax=ax4, colors='w', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False)
    ax4.quiver(data.lon[21:54:3], data.lat[10:55:1], u_after[3,10:55:1,21:54:3], v_after[3,10:55:1,21:54:3], headlength=5, headwidth=5)
    ax4.grid(True,linestyle=':')
    ax4.set_xlabel('Longitude')
    ax4.set_xlim(60,150)
    ax4.set_xticks(range(60,151,15))
    
    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0.1, hspace=0.1)
    #Colorbars
    cb1=f.colorbar(f1, ax=[ax1,ax2], use_gridspec=True, orientation = 'vertical',fraction=0.15, aspect=30)
    cb1.set_label('Precipitation, mm/day')
    cb2=f.colorbar(f3, ax=[ax3,ax4], use_gridspec=True, orientation = 'vertical',fraction=0.15, aspect=30)
    cb2.set_label('Surface Temperature, K')
    
    plt.savefig(plot_dir+'clim_before_after.pdf', format='pdf')
    plt.close()
    
    


clim_ba('ap_2', [30,39])
clim_ba('full_qflux', [18,39], land=True)
clim_ba('flat_qflux', [18,44], land=True)


