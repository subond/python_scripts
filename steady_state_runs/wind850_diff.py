"""
30/11/2017
Plot 850 hPa wind vectors and SST differences as a) heating strengthens and b) heating moves north
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr, model_constants as mc
from pylab import rcParams
import os

ctrl_90 = False

def get_wind(run1, run2):
        
    #Load data
    data_1 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run1+'.nc')
    data_2 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run2+'.nc')
    
    u = data_1.ucomp.sel(pfull=850)
    v = data_1.vcomp.sel(pfull=850)
    sst = data_1.t_surf
    u_diff = data_2.ucomp.sel(pfull=850) - u
    v_diff = data_2.vcomp.sel(pfull=850) - v
    sst_diff  = data_2.t_surf - sst
        
    return u_diff, v_diff, sst_diff, u, v, sst


def plot_diff(diff_type = 'strength'):
    
    rcParams['figure.figsize'] = 9, 12
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/steady_state_runs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    levels = np.arange(-2.,2.1,0.25)
    
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), 
         (ax10, ax11, ax12), (ax13, ax14, ax15), (ax16, ax17, ax18)) = plt.subplots(6, 3, sharex='col', sharey='row')
               
    plt.set_cmap('RdBu_r')

    ax_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17, ax18]
    lat_list = [0, 5, 10, 15, 20, 25]
    mag_list = [100, 200, 300]
    
    
    for i in range(18):
        
        mag = mag_list[np.mod(i,3)]
        lat = lat_list[i//3]
        #print mag, lat
        ax = ax_list[i]
        
        if diff_type == 'strength':
            run_1 = 'qflux_' + str(lat) + '_100'
            #print run_1
            run_2 = 'qflux_' + str(lat) + '_' + str(mag)
            #print run_2
        
        elif diff_type == 'loc':
            run_1 = 'qflux_' + '0_' + str(mag)
            #print run_1
            run_2 = 'qflux_' + str(lat) + '_' + str(mag)
            #print run_2
        
        else:
            raise ValueError('diff_type must be strength or loc')
        
        u_diff, v_diff, sst_diff, u, v, sst = get_wind(run_1, run_2)
        
        
        qflux_name = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/monsoon_ss_conts/input/qflux_95_lat_' + str(lat) + '_a_' + str(mag) + '.nc'
        data_qflux = xr.open_dataset(qflux_name)
        
        if (diff_type =='strength' and i in range(0,18,3)) or (diff_type =='loc' and i in range(3)):
            f1 = sst.plot.contourf(x='lon', y='lat', ax=ax, extend = 'both', add_labels=False, add_colorbar=False, levels=np.arange(270.,311.,5.))
            b = ax.quiver(u.lon[::5], v.lat[::2], u[::2,::5], v[::2,::5], scale=200., angles='xy')#,headwidth=5)
        else:
            f2 = sst_diff.plot.contourf(x='lon', y='lat', levels=np.arange(-10.,11.,1.), ax=ax, extend = 'both', add_labels=False, add_colorbar=False)
            b = ax.quiver(u_diff.lon[::5], v_diff.lat[::2], u_diff[::2,::5], v_diff[::2,::5], scale=200., angles='xy')#,headwidth=5)
            
        data_qflux.ocean_qflux.plot.contour(x='lon', y='lat', ax=ax, add_labels=False, colors='k', levels = np.arange(50.,350.,100.))
        ax.add_patch(patches.Rectangle((0,-60), 60, 20, facecolor="white"))
        b = ax.quiverkey(b, 0.075,0.05,10,'10 m/s', fontproperties={'weight': 'bold', 'size': 6}, color='k', labelcolor='k', labelsep=0.03)
        ax.set_yticks(range(-60,61,30))
        ax.set_ylim(-60.,60.)
        ax.set_xticks(range(0,361,90))
        ax.grid(True,linestyle=':')
        
    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.05, hspace=0.25, wspace=0.15)
    
    if diff_type == 'strength':
        ax_sublist = [ax2,ax3,ax5,ax6,ax8,ax9,ax11,ax12,ax14,ax15,ax17,ax18]
        cb1=f.colorbar(f2, ax=(ax_sublist), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)
        cb2=f.colorbar(f1, ax=([ax1,ax4,ax7,ax10,ax13,ax16]), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30)
        for label in cb2.ax.xaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
    else:
        ax_sublist = ax_list[3:]
        cb1=f.colorbar(f2, ax=(ax_sublist), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)
        cb2=f.colorbar(f1, ax=([ax1,ax2,ax3]), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30)
        for label in cb2.ax.xaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
            
    figname = 'uv_850_sst_diff_' + diff_type + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()


plot_diff(diff_type='loc')