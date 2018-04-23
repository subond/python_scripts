"""
Evaluate and plot momentum budget at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams

def abs_vort_hm(run, lev=150, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 10, 15
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    
    #Coriolis
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    abs_vort = (f + data.vor.sel(pfull=lev))*86400.
    
    abs_vort_zmean = abs_vort.sel(lon=lons).mean('lon')
    
    abs_vort_min = abs_vort
    abs_vort_min[0:32,:,:] = -1.*abs_vort_min[0:32,:,:]
    abs_vort_min = abs_vort_min.sel(lon=lons).min('lon')
    abs_vort_min[0:32,:] = -1.*abs_vort_min[0:32,:]
    
    
    levels = np.arange(-20.,20.1,2.)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    # Two subplots
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    plt.set_cmap('RdBu_r')
    #First plot
    f1=abs_vort_zmean.plot.contourf(ax=ax1, x='xofyear', y='lat', extend='both', levels=levels, add_colorbar=False, add_labels=False)
    abs_vort_zmean.plot.contour(ax=ax1, x='xofyear', y='lat', extend='both', levels=[-1000.,0.,1000.], add_colorbar=False, add_labels=False)
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)

    ax1.grid(True,linestyle=':')
    
    #Second plot
    abs_vort_min.plot.contourf(ax=ax2, x='xofyear', y='lat', extend='both', levels=levels, add_colorbar=False, add_labels=False)
    abs_vort_min.plot.contour(ax=ax2, x='xofyear', y='lat', extend='both', levels=[-1000.,0.,1000.], add_colorbar=False, add_labels=False)
    ax2.grid(True,linestyle=':')
    ax2.set_xticks(tickspace)
    ax2.set_xticklabels(labels,rotation=25)
    ax1.set_ylabel('Latitude')
    ax2.set_ylim(-60,60)

    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
    #Colorbar
    cb1=fig.colorbar(f1, ax=[ax1,ax2], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
    cb1.set_label('$day^{-1}$')
    
    if lonin == [-1.,361.]:
        figname = 'abs_vort.pdf'
    else:
        figname = 'abs_vort_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
        
        
abs_vort_hm('ap_2')
abs_vort_hm('full_qflux')
abs_vort_hm('full_qflux', lonin=[60.,150.])


