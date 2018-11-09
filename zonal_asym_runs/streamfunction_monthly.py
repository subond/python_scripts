''' 
14/08/2018 Plot upper and lower level streamfunction for each month
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from windspharm.xarray import VectorWind

plot_dir = '/scratch/rg419/plots/zonal_asym_runs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
def horiz_streamfun_monthly(run, land_mask=None):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    data.coords['month'] = (data.xofyear - 1) //6 + 1 
    data = data.groupby('month').mean(('xofyear'))

    # Create a VectorWind instance to handle the computation
    w = VectorWind(data.ucomp.sel(pfull=150.), data.vcomp.sel(pfull=150.))
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi_150, vchi_150, upsi_150, vpsi_150 = w.helmholtz()
    vel_pot_150 = (vel_pot - vel_pot.mean('lon'))/10.**6
    streamfun_150 = (streamfun - streamfun.mean('lon'))/10.**6
    
    w = VectorWind(data.ucomp.sel(pfull=850.), data.vcomp.sel(pfull=850.))
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi_850, vchi_850, upsi_850, vpsi_850 = w.helmholtz()
    vel_pot_850 = (vel_pot - vel_pot.mean('lon'))/10.**6
    streamfun_850 = (streamfun - streamfun.mean('lon'))/10.**6
    
    # Set figure parameters
    rcParams['figure.figsize'] = 15, 8
    rcParams['font.size'] = 14
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        f1 = streamfun_150[i,:,:].plot.contourf(ax=ax, x='lon', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-30.,30.1,5.), extend='both')
        streamfun_850[i,:,:].plot.contour(ax=ax, x='lon', y='lat', add_labels=False, add_colorbar=False, colors='k', levels=np.arange(0,12.1,2.))
        ax.contour(streamfun_850.lon, streamfun_850.lat, streamfun_850[i,:,:], colors='k', ls='--', levels=np.arange(-12.,0.,2.))
        if not land_mask==None:
            land = xr.open_dataset(land_mask)
            land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', alpha=0.5)
        #streamfun_850[i,:,:].plot.contour(ax=ax, x='lon', y='lat', add_labels=False, add_colorbar=False, cmap='PRGn', levels=np.arange(-12.,12.1,2.))

        i=i+1
        ax.set_ylim(-35,35)
        ax.set_yticks(np.arange(-30,31,15))
        ax.set_xlim(0,360)
        ax.set_xticks(np.arange(0,361,60))
        ax.grid(True,linestyle=':')
    
    for ax in [ax1,ax5,ax9]:
        ax.set_ylabel('Latitude')
    
    for ax in [ax9,ax10,ax11,ax12]:
        ax.set_xlabel('Longitude')
    
    plt.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Horizontal streamfunction')
    
    plt.savefig(plot_dir + 'streamfun_' + run + '.pdf', format='pdf')
        
    plt.close()
    
    
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        f1 = vel_pot_150[i,:,:].plot.contourf(ax=ax, x='lon', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-30.,30.1,5.), extend='both')
        vel_pot_850[i,:,:].plot.contour(ax=ax, x='lon', y='lat', add_labels=False, add_colorbar=False, colors='k', levels=np.arange(0,12.1,2.))
        ax.contour(vel_pot_850.lon, vel_pot_850.lat, vel_pot_850[i,:,:], colors='k', ls='--', levels=np.arange(-12.,0.,2.))
        if not land_mask==None:
            land = xr.open_dataset(land_mask)
            land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', alpha=0.5)
        #vel_pot_850[i,:,:].plot.contour(ax=ax, x='lon', y='lat', add_labels=False, add_colorbar=False, cmap='PRGn', levels=np.arange(-12.,12.1,2.))

        i=i+1
        ax.set_ylim(-35,35)
        ax.set_yticks(np.arange(-30,31,15))
        ax.set_xlim(0,360)
        ax.set_xticks(np.arange(0,361,60))
        ax.grid(True,linestyle=':')
    
    for ax in [ax1,ax5,ax9]:
        ax.set_ylabel('Latitude')
    
    for ax in [ax9,ax10,ax11,ax12]:
        ax.set_xlabel('Longitude')
    
    plt.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Velocity potential')
    
    plt.savefig(plot_dir + 'vel_pot_' + run + '.pdf', format='pdf')
        
    plt.close()
    
    
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        f1 = vel_pot_150[i,:,:].plot.contourf(ax=ax, x='lon', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-30.,30.1,5.), extend='both')
        b = axes[i].quiver(data.lon[::5], data.lat[::2], uchi_150[i,::2,::5], vchi_150[i,::2,::5], scale=100., angles='xy')
        if not land_mask==None:
            land = xr.open_dataset(land_mask)
            land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', alpha=0.5)
        i=i+1
        ax.set_ylim(-35,35)
        ax.set_yticks(np.arange(-30,31,15))
        ax.set_xlim(0,360)
        ax.set_xticks(np.arange(0,361,60))
        ax.grid(True,linestyle=':')
    
    for ax in [ax1,ax5,ax9]:
        ax.set_ylabel('Latitude')
    
    for ax in [ax9,ax10,ax11,ax12]:
        ax.set_xlabel('Longitude')
    
    plt.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Velocity potential')
    
    plt.savefig(plot_dir + 'vel_pot_vchi_' + run + '.pdf', format='pdf')
        
    plt.close()


horiz_streamfun_monthly('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
#horiz_streamfun_monthly('half_shallow_5', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
#horiz_streamfun_monthly('half_shallow_10', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
#horiz_streamfun_monthly('3q_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/3q_shallow.nc')
#horiz_streamfun_monthly('q_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/q_shallow.nc')
