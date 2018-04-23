"""
Dissect the components of the mean state terms in the vorticity budget
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams

def vort_eq_hm(run, lev=150, lonin=[-1.,361.]):
    
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/paper_1_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
        
    #Also load climatological data so that transient eddies can be calculated 
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_uv'+run+'.nc')
    print 'climatology loaded'
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
        
    dvordx = gr.ddx(vor)
    dvordy = gr.ddy(vor, vector=False)
        
    div = gr.ddx(data.ucomp) + gr.ddy(data.vcomp)
        
            
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    print 'starting plotting'
    
    # Easier to produce individual plots as scales are different.
    
    vor.sel(lon=lons).mean('lon').plot.contourf(x='xofyear', y='lat', extend = 'both', add_labels=False, levels=np.arange(-0.00024,0.00025,0.00002))
    plt.ylabel('Latitude')
    plt.ylim(-60,60)
    plt.yticks(np.arange(-60,61,30))
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.grid(True,linestyle=':')
    if lonin == [-1.,361.]:
        figname = 'vor_hm_' + run + '.pdf'
    else:
        figname = 'vor_hm_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    
    div.sel(lon=lons).mean('lon').plot.contourf(x='xofyear', y='lat', extend = 'both', add_labels=False, levels=np.arange(-0.000012,0.000012,0.000001))
    plt.ylabel('Latitude')
    plt.ylim(-60,60)
    plt.yticks(np.arange(-60,61,30))
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.grid(True,linestyle=':')
    if lonin == [-1.,361.]:
        figname = 'div_hm_' + run + '.pdf'
    else:
        figname = 'div_hm_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()


    dvordx.sel(lon=lons).mean('lon').plot.contourf(x='xofyear', y='lat', extend = 'both', add_labels=False)
    plt.ylabel('Latitude')
    plt.ylim(-60,60)
    plt.yticks(np.arange(-60,61,30))
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.grid(True,linestyle=':')
    if lonin == [-1.,361.]:
        figname = 'dvordx_hm_' + run + '.pdf'
    else:
        figname = 'dvordx_hm_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    
    dvordy.sel(lon=lons).mean('lon').plot.contourf(x='xofyear', y='lat', extend = 'both', add_labels=False, levels=np.arange(-1.2e-10,1.25e-10,0.1e-10))
    plt.ylabel('Latitude')
    plt.ylim(-60,60)
    plt.yticks(np.arange(-60,61,30))
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.grid(True,linestyle=':')
    if lonin == [-1.,361.]:
        figname = 'dvordy_hm_' + run + '.pdf'
    else:
        figname = 'dvordy_hm_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    
    data.ucomp.sel(lon=lons).mean('lon').plot.contourf(x='xofyear', y='lat', extend = 'both', add_labels=False)
    plt.ylabel('Latitude')
    plt.ylim(-60,60)
    plt.yticks(np.arange(-60,61,30))
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.grid(True,linestyle=':')
    if lonin == [-1.,361.]:
        figname = 'u_hm_' + run + '.pdf'
    else:
        figname = 'u_hm_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    
    data.vcomp.sel(lon=lons).mean('lon').plot.contourf(x='xofyear', y='lat', extend = 'both', add_labels=False, levels=np.arange(-15.,15.5,1.))
    plt.ylabel('Latitude')
    plt.ylim(-60,60)
    plt.yticks(np.arange(-60,61,30))
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.grid(True,linestyle=':')
    if lonin == [-1.,361.]:
        figname = 'v_hm_' + run + '.pdf'
    else:
        figname = 'v_hm_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    
      
vort_eq_hm('ap_2')
vort_eq_hm('full_qflux')
vort_eq_hm('full_qflux', lonin=[60.,150.])


