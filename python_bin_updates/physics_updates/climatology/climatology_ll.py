import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh


def plot_climatology(run, wind_and_rain=True, jets_and_temp=True, latent_and_lw=True, slp=True, land_mask=None):
    
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
    
    try:
        data['precipitation'] = (data.precipitation*86400.)
    except:
        data['precipitation'] = ((data.convection_rain + data.condensation_rain)*86400.)
    
    if wind_and_rain:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.precipitation[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(2.,15.,2.), add_labels=False, add_colorbar=False, extend='max', cmap='Blues')
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
        
        plt.savefig(plot_dir + 'wind_and_rain_' + run + '.pdf', format='pdf')
        plt.close()



    if jets_and_temp:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.t_surf[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels=np.arange(240.,300.,5.), add_labels=False, add_colorbar=False, extend='both')
            b = axes[i].quiver(data.lon[::5], data.lat[::2], data.ucomp.sel(pfull=150.,season=i)[::2,::5], data.vcomp.sel(pfull=150.,season=i)[::2,::5], scale=1000., angles='xy')
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
        
        plt.savefig(plot_dir + 'jets_and_temp_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    
    if latent_and_lw:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.flux_lhe[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(-300.,305.,20.), add_labels=False, add_colorbar=False, extend='both', cmap='RdBu')
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
        
        plt.savefig(plot_dir + 'latent_and_lw_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    if slp:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data.slp[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(970.,1021.,5.), add_labels=False, add_colorbar=False, extend='both')
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
        
        plt.savefig(plot_dir + 'slp_' + run + '.pdf', format='pdf')
        plt.close()
        
        

if __name__ == "__main__":
    
    plot_climatology('half_dry', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    plot_climatology('half_bright', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    #plot_climatology('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    #plot_climatology('q_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/q_shallow.nc')
    #plot_climatology('3q_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/3q_shallow.nc')
    