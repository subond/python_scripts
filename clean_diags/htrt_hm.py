"""
Evaluate and plot heat xrates at 850 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams
    
    
def mse_adv_hm(run, lev=850., filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 15, 10
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

    cp = 287.04/2.*7.
    L = 2.500e6
    data['mse_convection'] = (cp*data.dt_tg_convection + L*data.dt_qg_convection).sel(lon=lons).mean('lon')
    data['mse_condensation'] = (cp*data.dt_tg_condensation + L*data.dt_qg_condensation).sel(lon=lons).mean('lon')
    data['mse_diffusion'] = (cp*data.dt_tg_diffusion + L*data.dt_qg_diffusion).sel(lon=lons).mean('lon')
    data['mse_radiation'] = cp*data.tdt_rad.sel(lon=lons).mean('lon')

    data['mse'] = (cp*data.temp + L*data.sphum).sel(lon=lons).mean('lon')
    mse_max_loc = [data.lat[i] for i in data.mse.sel(pfull=lev).argmax('lat').values]

    levels = np.arange(-0.16,0.161,0.02)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    
    # Four subplots
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f1 = data.mse_convection.sel(pfull=lev).plot.contourf(x='xofyear', y='lat', ax=ax1, extend = 'both', add_colorbar=False, add_labels=False, cmap='RdBu_r', levels=levels)
    ax1.plot(data.xofyear,mse_max_loc, 'k')
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)
    ax1.grid(True,linestyle=':')
    
    #Second plot
    f2 = data.mse_condensation.sel(pfull=lev).plot.contourf(x='xofyear', y='lat', ax=ax2, extend = 'both', add_colorbar=False, add_labels=False, cmap='RdBu_r', levels=levels)
    ax2.plot(data.xofyear,mse_max_loc, 'k')
    ax2.grid(True,linestyle=':')
    
    #Third plot
    f3 = data.mse_diffusion.sel(pfull=lev).plot.contourf(x='xofyear', y='lat', ax=ax3, extend = 'both',  add_colorbar=False, add_labels=False, cmap='RdBu_r', levels=levels)
    ax3.plot(data.xofyear,mse_max_loc, 'k')
    ax3.grid(True,linestyle=':')
    ax3.set_ylabel('Latitude')
    ax3.set_xticks(tickspace)
    ax3.set_xticklabels(labels,rotation=25)    
    ax3.set_ylim(-60,60)
    
    #Fourth plot
    f4 = data.mse_radiation.sel(pfull=lev).plot.contourf(x='xofyear', y='lat', ax=ax4, extend = 'both', add_colorbar=False, add_labels=False, cmap='RdBu_r', levels=levels)
    ax4.plot(data.xofyear,mse_max_loc, 'k')
    ax4.grid(True,linestyle=':')
    ax4.set_xticks(tickspace)
    ax4.set_xticklabels(labels,rotation=25)
    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0.1, hspace=0.1)
    #Colorbar
    cb1=f.colorbar(f1, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30)
    cb1.set_label('$Jkg^{-1}s^{-1}$')
    
    if lonin == [-1.,361.]:
        figname = 'htrt_hm.pdf'
    else:
        figname = 'htrt_hm_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
        

mse_adv_hm('ap_2')
mse_adv_hm('full_qflux')
mse_adv_hm('full_qflux', lonin=[60.,150.])


