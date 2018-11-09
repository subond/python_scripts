'''11/10/2018 As gill_development.py but for data from JRA-55 and CMAP precip
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from windspharm.xarray import VectorWind


def plot_gill_dev(lev=85000, qscale=100., windtype='', ref_arrow=5):
    
    plot_dir = '/scratch/rg419/plots/zonal_asym_runs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data_u = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/ucomp_monthly/atmos_monthly_together.nc')
    data_v = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_monthly/atmos_monthly_together.nc')
    data_slp = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/slp_monthly/atmos_monthly_together.nc')
    data_precip = xr.open_dataset('/disca/share/rg419/CMAP_precip.mon.mean.nc')
    
    data_u = data_u['var33'].load().loc['1979-01':'2016-12'].groupby('time.month').mean('time')
    data_v = data_v['var34'].load().loc['1979-01':'2016-12'].groupby('time.month').mean('time')
    data_slp = data_slp['var2'].load().loc['1979-01':'2016-12'].groupby('time.month').mean('time') /100.
    data_precip = data_precip['precip'].load().loc['1979-01':'2016-12'].groupby('time.month').mean('time')
    
    data_u = data_u - data_u.mean('lon')
    data_v = data_v - data_v.mean('lon')
    data_slp = data_slp - data_slp.mean('lon')
    
    # Get rotational and divergent components of the flow
    w = VectorWind(data_u.sel(lev=lev), data_v.sel(lev=lev))
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    #print(uchi.lat)
    #print(data_u.lat)
    uchi_zanom = (uchi - uchi.mean('lon')).sortby('lat')
    vchi_zanom = (vchi - vchi.mean('lon')).sortby('lat')
    upsi_zanom = (upsi - upsi.mean('lon')).sortby('lat')
    vpsi_zanom = (vpsi - vpsi.mean('lon')).sortby('lat')
    
    
    # Start figure with 12 subplots
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 14
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6]
    titles = ['March','April', 'May', 'June', 'July', 'August']
    
    blist=[]
    for i in range(6):            
        f1 = data_precip[i+2,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(2.,15.,2.), add_labels=False, add_colorbar=False, extend='both',  cmap='Blues', zorder=1)
        #data_slp[i+3,:,:].plot.contour(x='lon', y='lat', ax=axes[i], levels = np.arange(-15.,16.,3.), add_labels=False, add_colorbar=False, extend='both', colors='0.5', alpha=0.5)
        axes[i].contour(data_slp.lon, data_slp.lat, data_slp[i+2,:,:], levels = np.arange(0.,16.,3.), colors='0.4', alpha=0.5, zorder=2)
        axes[i].contour(data_slp.lon, data_slp.lat, data_slp[i+2,:,:], levels = np.arange(-15.,0.,3.), colors='0.4', alpha=0.5, linestyle='--', zorder=2)
        if windtype=='div':
            b = axes[i].quiver(uchi_zanom.lon[::6], uchi_zanom.lat[::3], uchi_zanom[i+2,::3,::6], vchi_zanom[i+2,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
        elif windtype=='rot':
            b = axes[i].quiver(upsi_zanom.lon[::6], upsi_zanom.lat[::3], upsi_zanom[i+2,::3,::6], vpsi_zanom[i+2,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
        else:
            b = axes[i].quiver(data_u.lon[::6], data_u.lat[::3], data_u.sel(lev=lev)[i+2,::3,::6], data_v.sel(lev=lev)[i+2,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
        blist.append(b)
        axes[i].grid(True,linestyle=':')
        axes[i].set_ylim(-60.,60.)
        axes[i].set_yticks(np.arange(-60.,61.,30.))
        axes[i].set_xticks(np.arange(0.,361.,90.))
        axes[i].set_title(titles[i])
        #land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
        #land = xr.open_dataset(land_mask)
        #land.lsm[0,:,:].plot.contour(ax=axes[i], x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    
    for ax in [ax1, ax4]:
        ax.set_ylabel('Latitude')
    for ax in [ax4, ax5, ax6]:
        ax.set_xlabel('Longitude')
    
    ax4.quiverkey(blist[3], 0.,-0.5, ref_arrow, str(ref_arrow) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.93, bottom=0.05, hspace=0.25, wspace=0.2)
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.15, aspect=60, shrink=0.5)
    
    levstr=''
    windtypestr=''
    if lev != 85000:
        levstr = '_' + str(lev)
    if windtype != '':
        windtypestr = '_' + windtype

    plt.savefig(plot_dir + 'wind_and_slp_zanom_obs' + levstr + windtypestr + '.pdf', format='pdf')
    plt.close()
        

if __name__ == "__main__":
    
    plot_gill_dev()
    plot_gill_dev(windtype='rot')
    plot_gill_dev(windtype='div')
    plot_gill_dev(lev=20000, qscale=150., ref_arrow=10)
    plot_gill_dev(lev=20000, windtype='rot', qscale=150., ref_arrow=10)
    plot_gill_dev(lev=20000, windtype='div', qscale=150., ref_arrow=10)
    
                            
                            