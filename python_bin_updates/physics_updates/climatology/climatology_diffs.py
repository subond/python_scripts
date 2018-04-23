# 5/01/2017 Plot differences between climatologies of 2 runs to identify effects of bits of continents

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from data_handling_updates import model_constants as mc, gradients as gr


def plot_climatology(run_1, run_2, plot_title, wind_and_rain=True, jets_and_temp=True, latent_and_lw=True, slp=True, land_mask_1=None, land_mask_2=None):
    
    rcParams['figure.figsize'] = 10, 8
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/climatology/diffs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data_1 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run_1 + '.nc')
    data_2 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run_2 + '.nc')
    
    # Take seasonal averages
    data_1.coords['season'] = np.mod(data_1.xofyear + 5., 72.) // 18. 
    data_2.coords['season'] = np.mod(data_2.xofyear + 5., 72.) // 18. 
    #print data.season
    data_1 = data_1.groupby('season').mean(('xofyear'))
    data_2 = data_2.groupby('season').mean(('xofyear'))
    
    try:
        data_1['precipitation'] = (data_1.precipitation*86400.)
    except:
        data_1['precipitation'] = ((data_1.convection_rain + data_1.condensation_rain)*86400.)
    
    try:
        data_2['precipitation'] = (data_2.precipitation*86400.)
    except:
        data_2['precipitation'] = ((data_2.convection_rain + data_2.condensation_rain)*86400.)
    
    data = data_1 - data_2
    
    if wind_and_rain:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.precipitation[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(-5.,5.1,0.5), add_labels=False, add_colorbar=False, extend='both', cmap='RdBu')
            b = axes[i].quiver(data.lon[::5], data.lat[::2], data.ucomp.sel(pfull=850.,season=i)[::2,::5], data.vcomp.sel(pfull=850.,season=i)[::2,::5], scale=100., angles='xy')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask_1==None and not land_mask_2==None:
                land_1 = xr.open_dataset(land_mask_1)
                land_2 = xr.open_dataset(land_mask_2)
                land_1.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='0.5')
                land_2.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'wind_and_rain_' + plot_title + '.pdf', format='pdf')
        plt.close()



    if jets_and_temp:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.t_surf[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels=np.arange(-10.,10.1,1.), add_labels=False, add_colorbar=False, extend='both')
            b = axes[i].quiver(data.lon[::5], data.lat[::2], data.ucomp.sel(pfull=150.,season=i)[::2,::5], data.vcomp.sel(pfull=150.,season=i)[::2,::5], scale=500., angles='xy')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask_1==None and not land_mask_2==None:
                land_1 = xr.open_dataset(land_mask_1)
                land_2 = xr.open_dataset(land_mask_2)
                land_1.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='0.5')
                land_2.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
                        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'jets_and_temp_' + plot_title + '.pdf', format='pdf')
        plt.close()
    
    
    
    if latent_and_lw:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.flux_lhe[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(-50.,51.,5.), add_labels=False, add_colorbar=False, extend='both', cmap='RdBu')
            data.flux_lw[i,:,:].plot.contour(x='lon', y='lat', ax=axes[i], levels = np.arange(-100.,105.,10.), add_labels=False, add_colorbar=False, extend='both', cmap='Greys')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask_1==None and not land_mask_2==None:
                land_1 = xr.open_dataset(land_mask_1)
                land_2 = xr.open_dataset(land_mask_2)
                land_1.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='0.5')
                land_2.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
                        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'latent_and_lw_' + plot_title + '.pdf', format='pdf')
        plt.close()
    
    
    
    if slp:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.slp[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(-15.,16.,1.), add_labels=False, add_colorbar=False, extend='both')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask_1==None and not land_mask_2==None:
                land_1 = xr.open_dataset(land_mask_1)
                land_2 = xr.open_dataset(land_mask_2)
                land_1.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='0.5')
                land_2.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'slp_' + plot_title + '.pdf', format='pdf')
        plt.close()
        
        

if __name__ == "__main__":
    
    
    #plot_climatology('idealised_2cont', 'idealised_1cont', '2cont_1cont_diff',
    #     land_mask_1 = '/scratch/rg419/Experiments/wavenumber_2/input/2hill.nc', land_mask_2 = '/scratch/rg419/Experiments/wavenumber_2/input/1hill.nc')
    
    #plot_climatology('control_qflux', 'no_americas', 'full_no_am_diff',
    #     land_mask_1 = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc', land_mask_2 = '/scratch/rg419/python_scripts/land_era/land_era_no_america.nc')

    plot_climatology('control_qflux', 'frozen_am_0.4', 'full_frozen_0.4_diff',
         land_mask_1 = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc', land_mask_2 = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    
    #plot_climatology('control_qflux', 'no_TIP', 'full_no_TIP_diff',
    #     land_mask_1 = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc', land_mask_2 = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')

         
    #plot_climatology('idealised_land', 'idealised_land_no_am', 'no_am_diff',
    #     land_mask_1 = '/scratch/rg419/Experiments/idealised_land/input/land_all.nc', land_mask_2 = '/scratch/rg419/Experiments/idealised_land/input/land_aus_asia_tibet.nc')
         
    #plot_climatology('idealised_land_no_am', 'idealised_land_no_aus', 'no_aus_diff',
    #     land_mask_1 = '/scratch/rg419/Experiments/idealised_land/input/land_aus_asia_tibet.nc', land_mask_2 = '/scratch/rg419/Experiments/idealised_land/input/land_africa_asia_tibet.nc')
         
    #plot_climatology('idealised_land_no_aus', 'idealised_land_no_af', 'no_af_diff',
    #     land_mask_1 = '/scratch/rg419/Experiments/idealised_land/input/land_africa_asia_tibet.nc', land_mask_2 = '/scratch/rg419/Experiments/idealised_land/input/land_asia_tibet.nc')
         
    #plot_climatology('idealised_land_no_af', 'idealised_land_no_sa', 'no_sa_diff',
    #     land_mask_1 = '/scratch/rg419/Experiments/idealised_land/input/land_asia_tibet.nc', land_mask_2 = '/scratch/rg419/Experiments/idealised_land/input/land_EA_tibet.nc')
    
    #plot_climatology('idealised_land_no_sa', 'idealised_land_ea_only', 'no_tibet_diff',
    #     land_mask_1 = '/scratch/rg419/Experiments/idealised_land/input/land_EA_tibet.nc', land_mask_2 = '/scratch/rg419/Experiments/idealised_land/input/land_EA.nc')
         
    #plot_climatology('idealised_land_no_af', 'idealised_land_in_only', 'no_asia_diff',
    #     land_mask_1 = '/scratch/rg419/Experiments/idealised_land/input/land_asia_tibet.nc', land_mask_2 = '/scratch/rg419/Experiments/idealised_land/input/land_india.nc')
