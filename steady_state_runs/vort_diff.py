"""
30/11/2017
Plot absolute vorticity at 150 hPa as a) heating strengthens and b) heating moves north
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

def get_abs_vort(run1, run2):
        
    #Load data
    data_1 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_'+run1+'.nc')
    data_2 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_'+run2+'.nc')
    
    def calc_vort(data):
        # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
        f = 2 * mc.omega * np.sin(data.lat *np.pi/180)
        v_dx = gr.ddx(data.vcomp)  # dvdx
        u_dy = gr.ddy(data.ucomp)  # dudy
        vor = v_dx - u_dy + f
        abs_vort = vor.sel(pfull=150)*86400.
        return abs_vort
    
    vort_1 = calc_vort(data_1)
    vort_2 = calc_vort(data_2)
    
    vort_diff = vort_2 - vort_1
    
    return vort_diff, vort_1
    

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
        
        vort_diff, vort = get_abs_vort(run_1, run_2)
        
        
        qflux_name = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/monsoon_ss_conts/input/qflux_95_lat_' + str(lat) + '_a_' + str(mag) + '.nc'
        data_qflux = xr.open_dataset(qflux_name)
        
        if (diff_type =='strength' and i in range(0,18,3)) or (diff_type =='loc' and i in range(3)):
            f1 = vort.plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,1.), ax=ax, extend = 'both', add_labels=False, add_colorbar=False)
        else:
            f2 = vort_diff.plot.contourf(x='lon', y='lat', levels=np.arange(-5.,5.2,0.5), ax=ax, extend = 'both', add_labels=False, add_colorbar=False)
    
        data_qflux.ocean_qflux.plot.contour(x='lon', y='lat', ax=ax, add_labels=False, colors='k', levels = np.arange(50.,350.,100.))
        ax.set_yticks(range(-60,61,30))
        ax.set_ylim(-60.,60.)
        ax.set_xticks(range(0,361,90))
        ax.grid(True,linestyle=':')
        
    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.05, hspace=0.25, wspace=0.15)
    
    if diff_type == 'strength':
        ax_sublist = [ax2,ax3,ax5,ax6,ax8,ax9,ax11,ax12,ax14,ax15,ax17,ax18]
        cb1=f.colorbar(f2, ax=(ax_sublist), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)
        cb2=f.colorbar(f1, ax=([ax1,ax4,ax7,ax10,ax13,ax16]), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30)
        for label in cb1.ax.xaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
    else:
        ax_sublist = ax_list[3:]
        cb1=f.colorbar(f2, ax=(ax_sublist), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)
        cb2=f.colorbar(f1, ax=([ax1,ax2,ax3]), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30)

    for label in cb2.ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
            
    figname = 'vort_diff_' + diff_type + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()


plot_diff(diff_type='loc')
plot_diff(diff_type='strength')