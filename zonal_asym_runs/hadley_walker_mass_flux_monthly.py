''' 
31/08/2018 Plot local mass vertical mass flux associated with divergent flow in Hadley and Walker cells
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import mass_streamfunction
from data_handling_updates import model_constants as mc, gradients as gr
from windspharm.xarray import VectorWind


def h_w_mass_flux_monthly(run, lev=500., dp=5000.):
    
    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_v_clim_alllevs.nc' )
    data_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_u_clim_alllevs.nc' )
        
    plot_dir = '/scratch/rg419/plots/overturning_monthly/era/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data.coords['month'] = (data.xofyear - 1) //6 + 1 
    data = data.groupby('month').mean(('xofyear'))
    
    # Create a VectorWind instance to handle the computation
    w = VectorWind(data.ucomp.sel(pfull=np.arange(50.,950.,50.)), data.vcomp.sel(pfull=np.arange(50.,950.,50.)))
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    
    coslat = np.cos(data.lat * np.pi/180)
    
    # Evaluate mass fluxes for the zonal and meridional components (Walker and Hadley) following Schwendike et al. 2014
    mass_flux_zon = (gr.ddx(uchi)).cumsum('pfull') * dp * coslat/ mc.grav
    mass_flux_merid = (gr.ddy(vchi)).cumsum('pfull') * dp * coslat/ mc.grav
    
    # Set figure parameters
    rcParams['figure.figsize'] = 15, 11
    rcParams['font.size'] = 14
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        f1 = mass_flux_merid.sel(pfull=lev)[i,:,:].plot.contourf(ax=ax, x='lon', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-0.0065,0.0066,0.001), extend='both')
        mass_flux_zon.sel(pfull=lev)[i,:,:].plot.contour(ax=ax, x='lon', y='lat', add_labels=False, colors='k', levels=np.arange(0.0005,0.0066,0.001))
        mass_flux_zon.sel(pfull=lev)[i,:,:].plot.contour(ax=ax, x='lon', y='lat', add_labels=False, colors='0.5', levels=np.arange(-0.0065,-0.00049,0.001))

        i=i+1
        ax.set_ylim(-60,60)
        ax.set_xticks(np.arange(0,361,90))
        ax.set_yticks(np.arange(-60,61,30))
        ax.grid(True,linestyle=':')
    
    
    plt.subplots_adjust(left=0.05, right=0.97, top=0.95, bottom=0.1, hspace=0.2, wspace=0.2)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Vertical mass flux associated with meridional circulation, kgm$^{-2}$s$^{-1}$')
    
    figname = plot_dir + 'h_w_' + run + '.pdf'
    plt.savefig(figname, format='pdf')
    plt.close()
    

#h_w_mass_flux_monthly('half_shallow')
#h_w_mass_flux_monthly('3q_shallow')
#h_w_mass_flux_monthly('half_shallow_5')
#h_w_mass_flux_monthly('half_shallow_10')
h_w_mass_flux_monthly('half_nh_shallow')
h_w_mass_flux_monthly('half_10_shallow')
h_w_mass_flux_monthly('half_30_shallow')
h_w_mass_flux_monthly('nh_shallow')
h_w_mass_flux_monthly('10_shallow')
h_w_mass_flux_monthly('30_shallow')


#h_w_mass_flux_monthly('q_shallow')

