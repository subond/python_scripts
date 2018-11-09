''' 
04/09/2018 Plot local mass vertical mass flux associated with divergent flow in Hadley and Walker cells
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import mass_streamfunction
from data_handling_updates import model_constants as mc, gradients as gr
from windspharm.xarray import VectorWind


def h_w_mass_flux_monthly(lev=500., dp=5000.):
    
    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_v_clim_alllevs.nc' )
    data_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_u_clim_alllevs.nc' )
        
    plot_dir = '/scratch/rg419/plots/overturning_monthly/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    month_lengths = [31,28,31,30,31,30,31,31,30,31,30,31]
    month = []
    for i in range(1,13):
        month = month + [i]*month_lengths[i-1]
        
    data.coords['month'] = data.day_of_yr*0. + np.array(month)
    data = data.groupby('month').mean('day_of_yr')   
    
    data_u.coords['month'] = data_u.day_of_yr*0. + np.array(month)
    data_u = data_u.groupby('month').mean('day_of_yr')
    
    # Create a VectorWind instance to handle the computation
    w = VectorWind(data_u.u.sel(pfull=np.arange(50.,950.,50.)), data.v.sel(pfull=np.arange(50.,950.,50.)))
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    
    coslat = np.cos(data.lat * np.pi/180)
    
    # Evaluate mass fluxes for the zonal and meridional components (Walker and Hadley) following Schwendike et al. 2014
    mass_flux_zon = (gr.ddx(uchi)).cumsum('pfull') * dp * coslat/ mc.grav
    mass_flux_merid = (gr.ddy(vchi)).cumsum('pfull') * dp * coslat/ mc.grav
    
    land_mask='/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
    land = xr.open_dataset(land_mask)
    
    # Set figure parameters
    rcParams['figure.figsize'] = 15, 11
    rcParams['font.size'] = 18
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        f1 = mass_flux_merid.sel(pfull=lev)[i,:,:].plot.contourf(ax=ax, x='lon', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-0.0065,0.0066,0.001), extend='both')
        land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', ax=ax, levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        ax.fill_between([0,360], [-26, -26], [-7, -7], alpha=0.2, color='k')
        ax.fill_between([0,360], [7, 7], [26, 26], alpha=0.2, color='k')
        #mass_flux_zon.sel(pfull=lev)[i,:,:].plot.contour(ax=ax, x='lon', y='lat', add_labels=False, colors='k', levels=np.arange(0.0005,0.0066,0.001))
        #mass_flux_zon.sel(pfull=lev)[i,:,:].plot.contour(ax=ax, x='lon', y='lat', add_labels=False, colors='0.5', levels=np.arange(-0.0065,-0.00049,0.001))
        i=i+1
        ax.set_ylim(-60,60)
        ax.set_xticks(np.arange(0,361,90))
        ax.set_yticks(np.arange(-60,61,30))
        ax.grid(True,linestyle=':')
    
    
    plt.subplots_adjust(left=0.05, right=0.97, top=0.95, bottom=0.1, hspace=0.2, wspace=0.2)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Vertical mass flux associated with meridional circulation, kgm$^{-2}$s$^{-1}$')
    
    figname = plot_dir + 'hadley_era.pdf'
    plt.savefig(figname, format='pdf')
    plt.close()
    
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        f1 = mass_flux_zon.sel(pfull=lev)[i,:,:].plot.contourf(ax=ax, x='lon', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-0.0065,0.0066,0.001), extend='both')
        land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', ax=ax, levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        ax.fill_between([0,360], [-26, -26], [-7, -7], alpha=0.2, color='k')
        ax.fill_between([0,360], [7, 7], [26, 26], alpha=0.2, color='k')
        i=i+1
        ax.set_ylim(-60,60)
        ax.set_xticks(np.arange(0,361,90))
        ax.set_yticks(np.arange(-60,61,30))
        ax.grid(True,linestyle=':')
    
    
    plt.subplots_adjust(left=0.05, right=0.97, top=0.95, bottom=0.1, hspace=0.2, wspace=0.2)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Vertical mass flux associated with zonal circulation, kgm$^{-2}$s$^{-1}$')
    
    figname = plot_dir + 'walker_era.pdf'
    plt.savefig(figname, format='pdf')
    plt.close()


h_w_mass_flux_monthly()
