"""
Evaluate and plot momentum budget at 150 hPa before and after onset

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams


def partition_uv(data, lons, lev=150):
    
    uv_trans = (data.ucomp_vcomp - data.ucomp * data.vcomp).sel(pfull=lev)
    
    uv_stat = data.ucomp.sel(pfull=lev) * data.vcomp.sel(pfull=lev) - data.ucomp.sel(lon=lons, pfull=lev).mean('lon') * data.vcomp.sel(lon=lons, pfull=lev).mean('lon') 
    
    uv_tzav = data.ucomp.sel(lon=lons, pfull=lev).mean('lon') * data.vcomp.sel(lon=lons, pfull=lev).mean('lon') 

    data['uv_trans'] = (('lat','lon'), uv_trans)	
    data['uv_stat'] = (('lat','lon'), uv_stat)	
    data['uv_tzav']  = (('lat'), uv_tzav )
        
    
def uv_latlon(run, t, lev=150, land=False,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 10, 9
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        
    if land:
        land_data = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')
    
    t_dates = data.xofyear[t[0]:t[0]+4]
    data_before = data.sel(xofyear=t_dates).mean('xofyear')
    t_dates = data.xofyear[30:48]
    data_after = data.sel(xofyear=t_dates).mean('xofyear')
    
    #advective terms
    partition_uv(data_before, lons=lons, lev=lev)
    partition_uv(data_after, lons=lons, lev=lev)
    
    levels = np.arange(-200,200.1,20.)
    
    # Two subplots
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    
    #First plot
    f1 = data_before.uv_trans.plot.contourf(ax=ax1, x='lon', y='lat', levels = levels, add_colorbar=False, add_labels=False, extend='both')
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)
    ax1.set_xlim(0,150)
    ax1.grid(True,linestyle=':')
    
    #Secomd plot
    f2 = data_before.uv_stat.plot.contourf(ax=ax2, x='lon', y='lat', levels = levels, add_colorbar=False, add_labels=False, extend='both')
    ax2.grid(True,linestyle=':')
    ax2.invert_yaxis()
    ax2.set_ylim(-60,60)
    ax2.set_xlim(0,150)
    ax2.set_ylabel('Latitude')
    
    #Third plot
    f3 = data_after.uv_trans.plot.contourf(ax=ax3, x='lon', y='lat', levels = levels, add_colorbar=False, add_labels=False, extend='both')
    ax3.set_ylabel('Latitude')
    ax3.set_ylim(-60,60)
    ax3.set_xlim(0,150)
    ax3.grid(True,linestyle=':')
    ax3.set_ylabel('Latitude')
    ax3.set_xlabel('Longitude')
    
    #Fourth plot
    f4 = data_after.uv_stat.plot.contourf(ax=ax4, x='lon', y='lat', levels = levels, add_colorbar=False, add_labels=False, extend='both')
    ax4.grid(True,linestyle=':')
    ax4.invert_yaxis()
    ax4.set_xlim(0,150)
    ax4.set_ylim(-60,60)
    ax4.set_xlabel('Longitude')
    
    if land:
        for ax in [ax1,ax2,ax3,ax4]:
            land_data.land_mask.plot.contour(x='lon', y='lat', ax=ax, colors='k', linewidths=2, levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False)
    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
    #Colorbar
    cb1=f.colorbar(f1, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
    cb1.set_label('$ms^{-1}day^{-1}$')
    
    figname = 'uv_before_after.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
        

uv_latlon('full_qflux', [18,39], land=True)



