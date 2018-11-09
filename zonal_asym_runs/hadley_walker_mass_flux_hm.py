''' 
04/10/2018 Plot local mass vertical mass flux associated with divergent flow in Hadley and Walker cells as hovmoller
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import mass_streamfunction
from data_handling_updates import model_constants as mc, gradients as gr
from windspharm.xarray import VectorWind


def h_w_mass_flux_hm(run, dp=5000., lonin=[170,190]):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
        
    plot_dir = '/scratch/rg419/plots/zonal_asym_runs/' + run + '/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    # Create a VectorWind instance to handle the computation
    w = VectorWind(data.ucomp.sel(pfull=np.arange(50.,950.,50.)), data.vcomp.sel(pfull=np.arange(50.,950.,50.)))
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    
    coslat = np.cos(data.lat * np.pi/180)
    
    # Evaluate mass fluxes for the zonal and meridional components (Walker and Hadley) following Schwendike et al. 2014
    mass_flux_zon = (gr.ddx(uchi)).cumsum('pfull') * dp * coslat/ mc.grav
    mass_flux_merid = (gr.ddy(vchi)).cumsum('pfull') * dp * coslat/ mc.grav
    
    lats = [uchi.lat[i] for i in range(len(uchi.lat)) if uchi.lat[i] >= 0 and uchi.lat[i] <= 30]

    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    coslat = np.cos(uchi.lat*np.pi/180.)
    mass_flux_zon_areaweight = mass_flux_zon * coslat
    mass_flux_zon_mean = (mass_flux_zon_areaweight.sel(lat=lats).sum('lat') / coslat.sel(lat=lats).sum('lat')).sel(pfull=500)
    
    mass_flux_merid_mean = mass_flux_merid.sel(lon=lons).mean('lon').sel(pfull=500.)
    
    # Set figure parameters
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 14
    
    fig, (ax1, ax2) = plt.subplots(1, 2)
    
    f1 = mass_flux_zon_mean.plot.contourf(ax=ax1, x='lon', y='xofyear', add_labels=False, add_colorbar=False, levels=np.arange(-0.0065,0.0066,0.001), extend='both')
    f1 = mass_flux_merid_mean.plot.contourf(ax=ax2, x='xofyear', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-0.0065,0.0066,0.001), extend='both')
    ax2.set_ylim(-60,60)
    ax2.set_yticks(np.arange(-60,61,30))
    ax2.set_xticks(np.arange(0,73,12))
    ax1.set_yticks(np.arange(0,73,12))
    ax1.set_xticks(np.arange(0,361,90))
    ax1.grid(True,linestyle=':')
    ax2.grid(True,linestyle=':')
    
    plt.subplots_adjust(left=0.05, right=0.97, top=0.95, bottom=0.1, hspace=0.2, wspace=0.2)
    
    cb1=fig.colorbar(f1, ax=[ax1,ax2], use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Vertical mass flux, kgm$^{-2}$s$^{-1}$')
    
    figname = plot_dir + 'h_w_hm_' + run + '.pdf'
    plt.savefig(figname, format='pdf')
    plt.close()
    

h_w_mass_flux_hm('half_shallow')
#h_w_mass_flux_hm('3q_shallow')
#h_w_mass_flux_hm('half_shallow_5')
#h_w_mass_flux_hm('half_shallow_10')
#h_w_mass_flux_hm('half_nh_shallow')
#h_w_mass_flux_hm('half_10_shallow')
#h_w_mass_flux_hm('half_30_shallow')
#h_w_mass_flux_hm('nh_shallow')
#h_w_mass_flux_hm('10_shallow')
#h_w_mass_flux_hm('30_shallow')


#h_w_mass_flux_hm('q_shallow')

