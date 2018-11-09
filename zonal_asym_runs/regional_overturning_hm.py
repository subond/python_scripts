''' 
15/08/2018 Plot hm of local overturning streamfunction 
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import mass_streamfunction
from data_handling_updates import model_constants as mc
from windspharm.xarray import VectorWind


def overturning_hm(run, regions=[[350,10], [80,100], [170,190], [260,280]]):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
        
    plot_dir = '/scratch/rg419/plots/overturning_monthly/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    # Create a VectorWind instance to handle the computation
    w = VectorWind(data.ucomp.sel(pfull=np.arange(50.,950.,50.)), data.vcomp.sel(pfull=np.arange(50.,950.,50.)))
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    
    ds_chi = xr.Dataset({'vcomp': (vchi)},
                     coords={'xofyear': ('xofyear', vchi.xofyear),
                             'pfull': ('pfull', vchi.pfull),
                               'lat': ('lat', vchi.lat),
                               'lon': ('lon', vchi.lon)})
    
    def get_lons(lonin, data):
        if lonin[1]>lonin[0]:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
        else:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        return lons
    
    # Set figure parameters
    rcParams['figure.figsize'] = 10, 7
    rcParams['font.size'] = 14
    # Start figure with 4 subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    axes = [ax1,ax2,ax3,ax4]
    
    i=0
    for ax in axes:
        lons = get_lons(regions[i], data)
        psi_chi = mass_streamfunction(ds_chi, lons=lons, dp_in=50.)
        psi_chi /= 1.e9
        i=i+1
        f1=psi_chi.sel(pfull=500).plot.contourf(ax=ax, x='xofyear', y='lat', add_labels=False, add_colorbar=False, levels=np.arange(-500.,501.,100.), extend='both')
        
    ax1.set_title('West coast')
    ax2.set_title('Land')
    ax3.set_title('East coast')
    ax4.set_title('Ocean')
    
    for ax in [ax1,ax2,ax3,ax4]:
        ax.grid(True,linestyle=':')
        ax.set_ylim(-60,60)
        ax.set_yticks(np.arange(-60.,61.,30.))
        ax.set_xticks([0,18,36,54,72])
    
    ax3.set_xlabel('Pentad')
    ax4.set_xlabel('Pentad')
    ax1.set_ylabel('Latitude')
    ax3.set_ylabel('Latitude')
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Overturning streamfunction')
    
    # Save as a pdf
    plt.savefig(plot_dir + 'regional_overturning_hm_' + run + '.pdf', format='pdf')
    plt.close()
    
    data.close()



overturning_hm('half_shallow')
#overturning_hm('half_shallow_5')
#overturning_hm('half_shallow_10')
#overturning_hm('q_shallow', regions=[[350,10], [35,45], [80,100], [215,235]])
#overturning_hm('3q_shallow', regions=[[350,10], [125,145], [260,280], [305,325]])