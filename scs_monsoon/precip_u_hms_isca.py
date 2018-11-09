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

def precip_u_hms_isca(run, lonin=[110.,120.], rotfac=1.):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    lons = [data.lon[i].values for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

    (data.precipitation*86400.).sel(lon=lons).mean('lon').plot.contourf(ax=ax1, x='xofyear', y='lat', levels = np.arange(2.,15.,2.), add_colorbar=False, add_labels=False, cmap='Blues')
    (data.precipitation*86400.).sel(lon=lons).mean('lon').plot.contour(ax=ax1, x='xofyear', y='lat', levels = np.arange(6.,1006.,1000.), colors='k', add_labels=False)
    
    u_850 = data.ucomp.sel(pfull=850.).sel(lon=lons).mean('lon')
    v_850 = data.vcomp.sel(pfull=850.).sel(lon=lons).mean('lon')
    
    u_200 = data.ucomp.sel(pfull=200.).sel(lon=lons).mean('lon')
    v_200 = data.vcomp.sel(pfull=200.).sel(lon=lons).mean('lon')
    
    u_850.plot.contourf(ax=ax2, x='xofyear', y='lat', cmap = 'Greys', levels=np.arange(0.,1001.,1000.), add_labels=False, add_colorbar=False, extend='neither')
    b = ax2.quiver(data.xofyear, data.lat, u_850.T, v_850.T, scale=400.)
    ax2.add_patch(patches.Rectangle((0,-15), 7.5, 10, facecolor='0.9'))
    ax2.quiverkey(b, 0.05,0.05,10,'10 m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03)
    
    omega = 7.2921150e-5 * rotfac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vort = ((gr.ddx(data.vcomp.sel(pfull=200.)) - gr.ddy(data.ucomp.sel(pfull=200.))).sel(lon=lons).mean('lon') + f) * 86400.
    #vort = ((- gr.ddy(data.ucomp.sel(pfull=200.))).sel(lon=lons).mean('lon') + f) * 86400.
    
    vort.plot.contourf(ax=ax3, x='xofyear', y='lat', levels=np.arange(-12.,12.,2.), add_labels=False, add_colorbar=False)
    b = ax3.quiver(data.xofyear, data.lat, u_200.T, v_200.T, scale=1000.)
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
    
    plt.savefig(plot_dir+'precip_u_hms_' + run + '.pdf', format='pdf')
    plt.close()

#precip_u_hms_isca('idealised_2cont', lonin=[70.,90.])
precip_u_hms_isca('control_qflux')
precip_u_hms_isca('no_americas')
precip_u_hms_isca('no_TIP')
precip_u_hms_isca('control_qflux_0.500', rotfac=0.5)
precip_u_hms_isca('control_qflux_0.750', rotfac=0.75)
precip_u_hms_isca('control_qflux_1.250', rotfac=1.25)
precip_u_hms_isca('control_qflux_1.500', rotfac=1.5)
precip_u_hms_isca('control_qflux_1.750', rotfac=1.75)
precip_u_hms_isca('control_qflux_2.000', rotfac=2.)

precip_u_hms_isca('rt_0.500', rotfac=0.5)
precip_u_hms_isca('rt_0.750', rotfac=0.75)
precip_u_hms_isca('sn_1.000')
precip_u_hms_isca('rt_1.250', rotfac=1.25)
precip_u_hms_isca('rt_1.500', rotfac=1.5)
precip_u_hms_isca('rt_1.750', rotfac=1.75)
precip_u_hms_isca('rt_2.000', rotfac=2.)