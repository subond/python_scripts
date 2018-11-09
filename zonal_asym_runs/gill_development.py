'''11/10/2018 Plot monthly development of the Gill pattern and with sea level pressure and precipitation
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from windspharm.xarray import VectorWind
import matplotlib.patches as patches


def plot_gill_dev(run, land_mask=None, lev=850, qscale=100., windtype='', ref_arrow=5):
    
    plot_dir = '/scratch/rg419/plots/zonal_asym_runs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    # Take monthly averages
    data.coords['month'] = (data.xofyear - 1) //6 + 1 
    data = data.groupby('month').mean(('xofyear'))
    data_zanom = data - data.mean('lon')
    
    # Get rotational and divergent components of the flow
    w = VectorWind(data.ucomp.sel(pfull=lev), data.vcomp.sel(pfull=lev))
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    uchi_zanom = (uchi - uchi.mean('lon')).sortby('lat')
    vchi_zanom = (vchi - vchi.mean('lon')).sortby('lat')
    upsi_zanom = (upsi - upsi.mean('lon')).sortby('lat')
    vpsi_zanom = (vpsi - vpsi.mean('lon')).sortby('lat')
        
    data['precipitation'] = (data.precipitation*86400.)
    
    # Start figure with 12 subplots
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 14
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6]
    titles = ['May', 'June', 'July', 'August', 'September', 'October']
    
    blist=[]
    for i in range(6):    
        f1 = data.precipitation[i+4,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(2.,21.,2.), add_labels=False, add_colorbar=False, extend='both',  cmap='Blues', zorder=1)
        #data_zanom.slp[i+4,:,:].plot.contour(x='lon', y='lat', ax=axes[i], levels = np.arange(-15.,16.,3.), add_labels=False, colors='0.5', alpha=0.5)
        axes[i].contour(data_zanom.lon, data_zanom.lat, data_zanom.slp[i+4,:,:], levels = np.arange(0.,16.,3.), colors='0.4', alpha=0.5, zorder=2)
        axes[i].contour(data_zanom.lon, data_zanom.lat, data_zanom.slp[i+4,:,:], levels = np.arange(-15.,0.,3.), colors='0.4', alpha=0.5, linestyle='--', zorder=2)
        if windtype=='div':
            b = axes[i].quiver(data.lon[::6], data.lat[::3], uchi_zanom[i+4,::3,::6], vchi_zanom[i+4,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
        elif windtype=='rot':
            b = axes[i].quiver(data.lon[::6], data.lat[::3], upsi_zanom[i+4,::3,::6], vpsi_zanom[i+4,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
        else:
            b = axes[i].quiver(data.lon[::6], data.lat[::3], data_zanom.ucomp.sel(pfull=lev)[i+4,::3,::6], data_zanom.vcomp.sel(pfull=lev)[i+4,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
        blist.append(b)
        axes[i].grid(True,linestyle=':')
        axes[i].set_ylim(-60.,60.)
        axes[i].set_yticks(np.arange(-60.,61.,30.))
        axes[i].set_xticks(np.arange(0.,361.,90.))
        axes[i].set_title(titles[i])
        if not land_mask==None:
            land = xr.open_dataset(land_mask)
            land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    
    for ax in [ax1, ax4]:
        ax.set_ylabel('Latitude')
    for ax in [ax4, ax5, ax6]:
        ax.set_xlabel('Longitude')
    
    ax4.quiverkey(blist[3], 0.,-0.5, ref_arrow, str(ref_arrow) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.93, bottom=0.05, hspace=0.25, wspace=0.2)
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.15, aspect=60, shrink=0.5)
    
    levstr=''
    windtypestr=''
    if lev != 850:
        levstr = '_' + str(lev)
    if windtype != '':
        windtypestr = '_' + windtype

    plt.savefig(plot_dir + 'wind_and_slp_zanom_' + run + levstr + windtypestr + '.pdf', format='pdf')
    plt.close()
        

if __name__ == "__main__":
    
    plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='rot', lev=200., qscale=150., ref_arrow=10)
    plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='div', lev=200., qscale=150., ref_arrow=10)
    plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='rot')
    plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='div')
    plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', lev=200., qscale=150., ref_arrow=10)
    plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    plot_gill_dev('half_shallow_5', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    plot_gill_dev('half_shallow_10', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
    plot_gill_dev('half_nh_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_shallow.nc')
    plot_gill_dev('q_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/q_shallow.nc')
    plot_gill_dev('3q_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/3q_shallow.nc')
    

                            
                            