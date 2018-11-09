# 11/12/2017 Plot parts of evaporative flux to try to see where structure comes from.

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh


def plot_climatology_zanoms(run, wind_and_rain=True, jets_and_temp=True, latent_and_lw=True, slp=True, land_mask=None):
    
    rcParams['figure.figsize'] = 10, 8
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/climatology/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    # Take seasonal averages
    data.coords['season'] = np.mod(data.xofyear + 5., 72.) // 18. 
    #print data.season
    data = data.groupby('season').mean(('xofyear'))
    data = data - data.mean('lon')
    
    try:
        data['precipitation'] = (data.precipitation*86400.)
    except:
        data['precipitation'] = ((data.convection_rain + data.condensation_rain)*86400.)
    
    if wind_and_rain:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.precipitation[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(-10.,11.,1.), add_labels=False, add_colorbar=False, extend='both')
            b = axes[i].quiver(data.lon[::5], data.lat[::2], data.ucomp.sel(pfull=850.,season=i)[::2,::5], data.vcomp.sel(pfull=850.,season=i)[::2,::5], scale=200., angles='xy')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask==None:
                land = xr.open_dataset(land_mask)
                land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'wind_and_rain_zanom_' + run + '.pdf', format='pdf')
        plt.close()



    if jets_and_temp:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.t_surf[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(-15., 16.,1.), add_labels=False, add_colorbar=False, extend='both')
            b = axes[i].quiver(data.lon[::5], data.lat[::2], data.ucomp.sel(pfull=150.,season=i)[::2,::5], data.vcomp.sel(pfull=150.,season=i)[::2,::5], scale=500., angles='xy')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask==None:
                land = xr.open_dataset(land_mask)
                land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'jets_and_temp_zanom_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    
    if latent_and_lw:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.flux_lhe[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(-200.,201.,20.),  add_labels=False, add_colorbar=False, extend='both', cmap='RdBu')
            data.flux_lw[i,:,:].plot.contour(x='lon', y='lat', ax=axes[i], levels = np.arange(-500.,505.,20.), add_labels=False, add_colorbar=False, extend='both', cmap='Greys')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask==None:
                land = xr.open_dataset(land_mask)
                land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'latent_and_lw_zanom_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    if slp:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.slp[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(-15.,16.,1.), add_labels=False, add_colorbar=False, extend='both')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask==None:
                land = xr.open_dataset(land_mask)
                land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'slp_zanom_' + run + '.pdf', format='pdf')
        plt.close()
        
        

if __name__ == "__main__":
    
    #plot_climatology('continent_parts_no_am', land_mask = '/scratch/rg419/python_scripts/land_era/land_era_no_america.nc')
    plot_climatology_zanoms('half_dry', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    plot_climatology_zanoms('half_bright', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    #plot_climatology_zanoms('q_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/q_shallow.nc')
    #plot_climatology_zanoms('3q_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/3q_shallow.nc')
    
    #plot_climatology_zanoms('full_qflux', land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc', slp=False)
    
    #plot_climatology('idealised_1cont', land_mask = '/scratch/rg419/Experiments/wavenumber_2/input/1hill.nc')
    #plot_climatology('idealised_2cont', land_mask = '/scratch/rg419/Experiments/wavenumber_2/input/2hill.nc')
    
    #plot_climatology('no_americas', land_mask = '/scratch/rg419/python_scripts/land_era/land_era_no_america.nc')
    #plot_climatology('frozen_am', land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    #plot_climatology('no_TIP', land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    #plot_climatology_zanoms('control_qflux', land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
                            
                            