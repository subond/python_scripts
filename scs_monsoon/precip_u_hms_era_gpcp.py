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

rcParams['figure.figsize'] = 12, 13
rcParams['font.size'] = 14

plot_dir = '/scratch/rg419/plots/scs_monsoon/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

def precip_u_hms_era(lonin=[110.,120.]):
    
    data_precip = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/gpcp_detrendedpentads.nc')
    data_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_850_clim.nc')
    data_v = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_v/era_v_850_clim.nc')
    data_vo = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_vo/era_vo_200_clim.nc')
    
    
    lons_precip = [data_precip.lon[i].values for i in range(len(data_precip.lon)) if data_precip.lon[i] >= lonin[0] and data_precip.lon[i] < lonin[1]]
    lons_era = [data_u.lon[i].values for i in range(len(data_u.lon)) if data_u.lon[i] >= lonin[0] and data_u.lon[i] < lonin[1]]

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

    data_precip.precip_clim.sel(lon=lons_precip).mean('lon').plot.contourf(ax=ax1, x='xofyear', y='lat', levels = np.arange(2.,15.,2.), add_colorbar=False, add_labels=False, cmap='Blues')
    data_precip.precip_clim.sel(lon=lons_precip).mean('lon').plot.contour(ax=ax1, x='xofyear', y='lat', levels = np.arange(6.,1006.,1000.), colors='k', add_labels=False)
    
    u_850 = data_u.u_850.sel(lon=lons_era).mean('lon')
    v_850 = data_v.v_850.sel(lon=lons_era).mean('lon')
    u_850.coords['xofyear'] = np.mod( u_850.day_of_yr, 365.) //5 + 1.  
    u_850 = u_850.groupby('xofyear').mean(('day_of_yr'))
    v_850.coords['xofyear'] = np.mod( v_850.day_of_yr, 365.) //5 + 1.  
    v_850 = v_850.groupby('xofyear').mean(('day_of_yr'))
        
    u_850.plot.contourf(ax=ax2, x='xofyear', y='lat', cmap = 'Greys', levels=np.arange(0.,1001.,1000.), add_labels=False, add_colorbar=False, extend='neither')    
    b = ax2.quiver(u_850.xofyear, u_850.lat, u_850, v_850, scale=400.)
    ax2.add_patch(patches.Rectangle((0,-15), 7.5, 10, facecolor='0.9'))
    ax2.quiverkey(b, 0.05,0.05,10,'10 m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03)
    
    omega = 7.2921150e-5 
    f = 2 * omega * np.sin(data_vo.lat *np.pi/180)
    vo_200 = (data_vo.vo_200.sel(lon=lons_era).mean('lon') + f) * 86400.
    vo_200.coords['xofyear'] = np.mod( vo_200.day_of_yr, 365.) //5 + 1.  
    vo_200 = vo_200.groupby('xofyear').mean(('day_of_yr'))
    
    data_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_200_clim.nc')
    data_v = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_v/era_v_200_clim.nc')
    u_200 = data_u.u_200.sel(lon=lons_era).mean('lon')
    v_200 = data_v.v_200.sel(lon=lons_era).mean('lon')
    u_200.coords['xofyear'] = np.mod( u_200.day_of_yr, 365.) //5 + 1.  
    u_200 = u_200.groupby('xofyear').mean(('day_of_yr'))
    v_200.coords['xofyear'] = np.mod( v_200.day_of_yr, 365.) //5 + 1.  
    v_200 = v_200.groupby('xofyear').mean(('day_of_yr'))
    
    vo_200.plot.contourf(ax=ax3, x='xofyear', y='lat', levels=np.arange(-12.,12.,2.), add_labels=False, add_colorbar=False)
    b = ax3.quiver(u_200.xofyear, u_200.lat, u_200, v_200, scale=1000.)
    ax3.add_patch(patches.Rectangle((0,-15), 7.5, 10, facecolor='0.9'))
    ax3.quiverkey(b, 0.05,0.05,50,'50 m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03)
    
    ax3.set_xlabel('Pentad')
    ax1.set_ylabel('Latitude')
    ax2.set_ylabel('Latitude')
    ax3.set_ylabel('Latitude')
    ax1.set_ylim([-15,45])
    ax2.set_ylim([-15,45])
    ax3.set_ylim([-15,45])
    ax1.grid(True,linestyle=':', linewidth=2, color='k')
    ax2.grid(True,linestyle=':', linewidth=2, color='k')
    ax3.grid(True,linestyle=':', linewidth=2, color='k')
    
    plt.savefig(plot_dir+'precip_u_hms_era_' + str(lonin[0]) + '_' + str(lonin[1]) + '.pdf', format='pdf')
    plt.close()


precip_u_hms_era(lonin=[60.,150.])