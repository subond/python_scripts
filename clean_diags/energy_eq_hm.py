"""
Load in the vorticity budget terms and produce a lat-time plot at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams

def energy_eq_hm(run, lev=150, lonin=[-1.,361.], levels = np.arange(-0.01,0.011,0.001)):
    
    rcParams['figure.figsize'] = 10, 9
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Load in energy budget term means
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/energy_eq_'+run+'.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    pe_to_ke_hm = data.pe_to_ke.sel(pfull=lev, lon=lons).mean('lon')
    conv_1_hm = data.conv_1.sel(pfull=lev, lon=lons).mean('lon')
    conv_2_hm = data.conv_2.sel(pfull=lev, lon=lons).mean('lon')
    conv_3_hm = data.conv_3.sel(pfull=lev, lon=lons).mean('lon')
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    # Four subplots
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    
    #First plot
    f2 = pe_to_ke_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax1, extend = 'both', add_labels=False, add_colorbar=False)
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)
    ax1.grid(True,linestyle=':')

    #Second plot
    conv_1_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax2, extend = 'both', add_labels=False, add_colorbar=False)
    ax2.set_ylim(-60,60)
    ax2.grid(True,linestyle=':')
    
    #Third plot
    conv_3_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax3, extend = 'both', add_labels=False, add_colorbar=False)
    ax3.set_ylabel('Latitude')
    ax3.set_ylim(-60,60)
    ax3.set_xticks(tickspace)
    ax3.set_xticklabels(labels,rotation=25)
    ax3.grid(True,linestyle=':')
    
    #Fourth plot
    conv_2_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax4, extend = 'both', add_labels=False, add_colorbar=False)
    ax4.set_ylim(-60,60)
    ax4.set_xticks(tickspace)
    ax4.set_xticklabels(labels,rotation=25)
    ax4.grid(True,linestyle=':')
    
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.3, wspace=0.1)
    #Colorbar
    cb1=f.colorbar(f2, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30, shrink=0.9)
    cb1.set_label('Energy conversion, J/s')
    
    
    if lonin == [-1.,361.]:
        figname = 'energy_budg_hm_' + str(int(lev)) + '.pdf'
    else:
        figname = 'energy_budg_hm_' + str(int(lev)) + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    
#energy_eq_hm('ap_2')
#energy_eq_hm('ap_2', lev=850., levels=np.arange(-0.006,0.007,0.001))
energy_eq_hm('full_qflux')
energy_eq_hm('full_qflux', lonin=[60.,150.])
energy_eq_hm('full_qflux', lev=850., levels=np.arange(-0.001,0.0011,0.0001))
energy_eq_hm('full_qflux', lonin=[60.,150.], lev=850., levels=np.arange(-0.001,0.0011,0.0001))

