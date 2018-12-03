""" 03/12/2018 Plot histograms of the onset dates in isca and in jra data and save
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sh
from pylab import rcParams
from climatology import scsm_onset

rcParams['figure.figsize'] = 12, 6
rcParams['font.size'] = 14

plot_dir = '/scratch/rg419/plots/scs_monsoon/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

data = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_daily/atmos_daily_together.nc', chunks={'time': 30})
data = data.sel(lev=85000.)
data.to_netcdf('/disca/share/rg419/jra_vcomp_daily_850.nc')


def pentad_mean_climatology(data, years):  # Function to get pentad of year
    pentad_years = np.array([])
    for year in years:
        data_year = data.sel(time=str(year))
        if len(data_year.time)==366:
            pentad = np.repeat(np.arange(1., 74.), 5)
            pentad = np.insert(pentad, 10, 2)    
        else:
            pentad = np.repeat(np.arange(1., 74.), 5)    
        np.concatenate((pentad_years, pentad))
        print(pentad_years.shape)
        
    data = data.assign_coords(pentad = ('time', pentad_years))
    print(data.pentad_years.values[0:720])
    
    data = data.groupby('pentad').mean(('time'))
    
    return data_pentads


def precip_u_hms_jra(lonin=[110.,120.]):
    
    data_precip = xr.open_dataset('/disca/share/rg419/CMAP_precip.pentad.mean.nc')
    #data_u = xr.open_dataset('/disca/share/rg419/jra_ucomp_daily_200.nc')
    data_u = xr.open_dataset('/disca/share/rg419/jra_ucomp_daily_850.nc')
    data_u = data_u['var33'].load().loc['1958-01':'2016-12']
    #data_u_clim = data_u.groupby('time.dayofyear').mean('time')
    #print(data_u_clim)

    # v has different time coord to u, presumably due to how Stephen has downloaded/averaged. I think the two are equivalent, so just substitute the time dimension into v
    data_v_temp = xr.open_dataset('/disca/share/rg419/jra_vcomp_daily_200.nc')
    data_v = xr.DataArray(data_v_temp.var34.values, coords=data_u.coords, dims=data_u.dims)
    
    #data_u = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/ucomp_daily/atmos_daily_together.nc')
    #data_v = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_daily/atmos_daily_together.nc')
    
    print('files opened')
    
    data_precip.coords['pentad'] = (('time'), np.tile(np.arange(1,74),38))
    data_precip = data_precip.groupby('pentad').mean('time')
    
    data_u = pentad_mean_climatology(data_u, np.arange(1979,2017))
    data_v = pentad_mean_climatology(data_v, np.arange(1979,2017))
    
    print('pentad means taken, plotting')
    
    lons_precip = [data_precip.lon[i].values for i in range(len(data_precip.lon)) if data_precip.lon[i] >= lonin[0] and data_precip.lon[i] < lonin[1]]
    lons_jra = [data_u.lon[i].values for i in range(len(data_u.lon)) if data_u.lon[i] >= lonin[0] and data_u.lon[i] < lonin[1]]

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2 , sharex='col', sharey='row')

    data_precip.precip.sel(lon=lons_precip).mean('lon').plot.contourf(ax=ax1, x='pentad', y='lat', levels = np.arange(2.,15.,2.), add_colorbar=False, add_labels=False, cmap='Blues')
    data_precip.precip.sel(lon=lons_precip).mean('lon').plot.contour(ax=ax1, x='pentad', y='lat', levels = np.arange(6.,1006.,1000.), colors='k', add_labels=False)
    
    u_850 = data_u.var33.sel(lon=lons_jra, lev=85000.).mean('lon')
    v_850 = data_v.var34.sel(lon=lons_jra, lev=85000.).mean('lon')
        
    u_850.plot.contourf(ax=ax3, x='xofyear', y='lat', cmap = 'Greys', levels=np.arange(0.,1001.,1000.), add_labels=False, add_colorbar=False, extend='neither')    
    b = ax3.quiver(u_850.xofyear[::2], u_850.lat, u_850[:,::2], v_850[:,::2], scale=200.)
    
    print('JRA and CMAP plotted, open Isca data')
    
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

    
    plt.savefig(plot_dir+'precip_u_hms_jra_isca.pdf', format='pdf')
    plt.close()


#precip_u_hms_jra()
#precip_u_hms_jra(lonin=[60.,150.])