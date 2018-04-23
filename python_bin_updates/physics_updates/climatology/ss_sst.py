# 1/12/2017 Plot steady state run SSTs and precip

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from physics import model_constants as mc, gradients as gr


g = 9.8
cp = 287.04/2*7
L = 2.500e6

def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons


def ss_sst_plot(runs, lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 6, 6
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/surface_fluxes/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + runs[0] + '.nc')
    lons = pick_lons(data, lonin)
    
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    
    for run in runs:
        data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
        
        try:
            precip_plot = (data.precipitation*86400.).sel(lon=lons).mean('lon')
        except:
            precip_plot = ((data.convection_rain + data.condensation_rain)*86400.).sel(lon=lons).mean('lon')
        
        data = data.sel(lon=lons).mean('lon')
        
        f1 = data.t_surf.plot(ax=ax1)
        precip_plot.plot(ax=ax2)
        
        ax1.set_ylabel('SST, K')
        ax2.set_ylabel('Precip, mm/day')
        ax1.set_xlabel(' ')
        ax2.set_xlabel('Latitude')
        ax2.set_xlim([-60,60])
        ax1.set_ylim([240,320])
    
    plt.legend(runs, loc = 'upper left', prop={'size': 10})

    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1, hspace=0.2)

    plt.savefig(plot_dir + 'ss_sst_plot.pdf', format='pdf')
    plt.close()



def evolving_sst_plot(run, pentads=[40,42,44,46,48,50], lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 6, 6
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/surface_fluxes/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    lons = pick_lons(data, lonin)
    
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    
    for i in pentads:
        data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
        
        try:
            precip_plot = (data.precipitation*86400.).sel(lon=lons).mean('lon')
        except:
            precip_plot = ((data.convection_rain + data.condensation_rain)*86400.).sel(lon=lons).mean('lon')
        
        data = data.sel(lon=lons).mean('lon')
        
        f1 = data.t_surf.sel(xofyear=i).plot(ax=ax1)
        precip_plot.sel(xofyear=i).plot(ax=ax2)
        
        ax1.set_ylabel('SST, K')
        ax2.set_ylabel('Precip, mm/day')
        ax1.set_xlabel(' ')
        ax2.set_xlabel('Latitude')
        ax2.set_xlim([-60,60])
        ax1.set_ylim([240,320])
        ax2.set_ylim([0,20])
        ax1.set_title(' ')
        ax2.set_title(' ')
        
    
    plt.legend([str(pentad) for pentad in pentads], loc = 'upper left', prop={'size': 10})

    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1, hspace=0.2)

    plt.savefig(plot_dir + 'sst_plot_' + run + '.pdf', format='pdf')
    plt.close()
    
    

if __name__ == "__main__":
    
    #runs = ['ss_90.000', 'ss_95.000', 'ss_100.000', 'ss_105.000', 'ss_110.000', 'ss_115.000']
    #ss_sst_plot(runs)
    
    evolving_sst_plot('sn_1.000')
    

