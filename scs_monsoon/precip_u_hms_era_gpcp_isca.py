""" 14/06/2018 Plot histograms of the onset dates in isca and in era data and save
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sh
from pylab import rcParams
from climatology import scsm_onset
from data_handling_updates import isca_load_and_reshape, gradients as gr

rcParams['figure.figsize'] = 12, 6
rcParams['font.size'] = 14

plot_dir = '/scratch/rg419/plots/scs_monsoon/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

def precip_u_hms_era(lonin=[110.,120.]):
    
    data_precip = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/gpcp_detrendedpentads.nc')
    data_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_850_clim.nc')
    data_v = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_v/era_v_850_clim.nc')
    
    lons_precip = [data_precip.lon[i].values for i in range(len(data_precip.lon)) if data_precip.lon[i] >= lonin[0] and data_precip.lon[i] < lonin[1]]
    lons_era = [data_u.lon[i].values for i in range(len(data_u.lon)) if data_u.lon[i] >= lonin[0] and data_u.lon[i] < lonin[1]]

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2 , sharex='col', sharey='row')

    data_precip.precip_clim.sel(lon=lons_precip).mean('lon').plot.contourf(ax=ax1, x='xofyear', y='lat', levels = np.arange(2.,15.,2.), add_colorbar=False, add_labels=False, cmap='Blues')
    data_precip.precip_clim.sel(lon=lons_precip).mean('lon').plot.contour(ax=ax1, x='xofyear', y='lat', levels = np.arange(6.,1006.,1000.), colors='k', add_labels=False)
    
    u_850 = data_u.u_850.sel(lon=lons_era).mean('lon')
    v_850 = data_v.v_850.sel(lon=lons_era).mean('lon')
    u_850.coords['xofyear'] = np.mod( u_850.day_of_yr, 365.) //5 + 1.  
    u_850 = u_850.groupby('xofyear').mean(('day_of_yr'))
    v_850.coords['xofyear'] = np.mod( v_850.day_of_yr, 365.) //5 + 1.  
    v_850 = v_850.groupby('xofyear').mean(('day_of_yr'))
        
    u_850.plot.contourf(ax=ax3, x='xofyear', y='lat', cmap = 'Greys', levels=np.arange(0.,1001.,1000.), add_labels=False, add_colorbar=False, extend='neither')    
    b = ax3.quiver(u_850.xofyear[::2], u_850.lat, u_850[:,::2], v_850[:,::2], scale=200.)
    
    data_isca = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/sn_1.000.nc')
    
    (data_isca.precipitation*86400.).mean('lon').plot.contourf(ax=ax2, x='xofyear', y='lat', levels = np.arange(2.,15.,2.), add_colorbar=False, add_labels=False, cmap='Blues')
    (data_isca.precipitation*86400.).mean('lon').plot.contour(ax=ax2, x='xofyear', y='lat', levels = np.arange(6.,1006.,1000.), colors='k', add_labels=False)
    
    u_850 = data_isca.ucomp.sel(pfull=850.).mean('lon')
    v_850 = data_isca.vcomp.sel(pfull=850.).mean('lon')
    
    u_850.plot.contourf(ax=ax4, x='xofyear', y='lat', cmap = 'Greys', levels=np.arange(0.,1001.,1000.), add_labels=False, add_colorbar=False, extend='neither')
    b1 = ax4.quiver(data_isca.xofyear[::2], data_isca.lat, u_850.T[:,::2], v_850.T[:,::2], scale=200.)
    
    ax1.grid(True,linestyle=':', linewidth=2, color='k')
    ax2.grid(True,linestyle=':', linewidth=2, color='k')
    ax3.grid(True,linestyle=':', linewidth=2, color='k')
    ax4.grid(True,linestyle=':', linewidth=2, color='k')
    
    ax3.add_patch(patches.Rectangle((0,-15), 12, 10, facecolor='0.9'))
    ax3.quiverkey(b, 0.07,0.05,10,'10 m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03)
    ax4.add_patch(patches.Rectangle((0,-15), 12, 10, facecolor='0.9'))
    ax4.quiverkey(b1, 0.07,0.05,10,'10 m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03)
    
    ax3.set_xlabel('Pentad')
    ax4.set_xlabel('Pentad')
    ax1.set_ylabel('Latitude')
    ax3.set_ylabel('Latitude')
    ax1.set_ylim([-15,45])
    ax2.set_ylim([-15,45])
    ax3.set_ylim([-15,45])
    ax4.set_ylim([-15,45])

    
    plt.savefig(plot_dir+'precip_u_hms_era_isca.pdf', format='pdf')
    plt.close()


precip_u_hms_era(lonin=[60.,150.])