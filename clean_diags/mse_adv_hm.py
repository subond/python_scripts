"""
Evaluate and plot moist static energy budget at 850 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams
    
    
def mse_adv_hm(run, lev=850., filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 7, 5
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

    cp = 287.04/2*7
    L = 2.500e6
    data['mse_v'] = (cp*data.vcomp_temp + L*data.sphum_v).sel(lon=lons).mean('lon')
    data['mse_v_mean'] = cp*data.vcomp.sel(lon=lons).mean('lon') * data.temp.sel(lon=lons).mean('lon') + L*data.sphum.sel(lon=lons).mean('lon') * data.vcomp.sel(lon=lons).mean('lon')
    data['mse'] = (cp*data.temp + L*data.sphum).sel(lon=lons).mean('lon')
    mse_max_loc = [data.lat[i] for i in data.mse.sel(pfull=lev).argmax('lat').values]

    levels = np.arange(-8000,8000.1,1000.)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    data.mse_v.sel(pfull=lev).plot.contourf(x='xofyear', y='lat', levels=np.arange(-10.e5,10.1e5,20.e4), extend='both', add_labels=False)
    data.mse.sel(pfull=lev).plot.contour(x='xofyear', y='lat', colors='k', levels=np.arange(200000.,400000.,20000.), add_colorbar=False, add_labels=False)
    plt.plot(data.xofyear,mse_max_loc)
    plt.grid(True,linestyle=':')
    plt.xticks(tickspace,labels,rotation=25)
    plt.tight_layout()  
    plt.ylim(-60,60)
    plt.ylabel('Latitude')
    if lonin == [-1.,361.]:
        figname = 'mse_v.pdf'
    else:
        figname = 'mse_v_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()   

    
    data.mse_v_mean.sel(pfull=lev).plot.contourf(x='xofyear', y='lat', levels=np.arange(-10.e5,10.1e5,20.e4), extend='both', add_labels=False)
    data.mse.sel(pfull=lev).plot.contour(x='xofyear', y='lat', colors='k', levels=np.arange(200000.,400000.,20000.), add_colorbar=False, add_labels=False)
    plt.plot(data.xofyear,mse_max_loc)
    plt.grid(True,linestyle=':')
    plt.xticks(tickspace,labels,rotation=25)
    plt.tight_layout()  
    plt.ylim(-60,60)
    plt.ylabel('Latitude')
    if lonin == [-1.,361.]:
        figname = 'mse_v_mean.pdf'
    else:
        figname = 'mse_v_mean_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()   
        
#mse_budg_hm('ap_20_qflux')
#mse_budg_hm('ap_20_qflux', lonin=[60.,150.])
#mse_budg_hm('am_qflux')
#mse_budg_hm('am_qflux', lonin=[60.,150.])
mse_adv_hm('ap_2')
mse_adv_hm('full_qflux')
mse_adv_hm('full_qflux', lonin=[60.,150.])
#mom_budg_hm('flat_qflux', [121,481])
#mom_budg_hm('flat_qflux', [121,481], lonin=[60.,150.])


