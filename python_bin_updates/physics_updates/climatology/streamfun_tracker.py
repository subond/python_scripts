# 11/01/2018 Evaluate streamfunction, plot for JJA with centre by month on top

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from physics import model_constants as mc, gradients as gr
from windspharm.xarray import VectorWind


def plot_sf_vp(run, land_mask=None):
    
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/climatology/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    uwnd = data.ucomp.sel(pfull=200.)
    vwnd = data.vcomp.sel(pfull=200.)
    # Create a VectorWind instance to handle the computation
    w = VectorWind(uwnd, vwnd)

    # Compute variables
    streamfun, vel_pot = w.sfvp()

    def mn_sn_av(da):
        #Take seasonal and monthly averages
        da = da/10.**6
        da.coords['season'] = np.mod(da.xofyear + 5., 72.) // 18. 
        da.coords['month'] = (da.xofyear-1) //6 + 1     
        da_sn = da.groupby('season').mean(('xofyear'))
        da_mn = da.groupby('month').mean(('xofyear'))
        return da_sn, da_mn
    
    streamfun_sn, streamfun_mn = mn_sn_av(streamfun)
    vel_pot_sn, vel_pot_mn = mn_sn_av(vel_pot)
        

    
    lats= [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= 5.]    
    lons= [data.lon[i] for i in range(len(data.lon)) if data.lon[i] <= 175. and data.lon[i] >=25.]    
    
    sf_max_lat = np.zeros((12,))
    vp_min_lat = np.zeros((12,))
    sf_max_lon = np.zeros((12,))
    vp_min_lon = np.zeros((12,))
    
    for i in range(12):

        streamfun_asia_i = streamfun_mn[i,:,:].sel(lat=lats, lon=lons)
        sf_max_i = streamfun_asia_i.where(streamfun_asia_i==streamfun_asia_i.max(), drop=True)
        
        sf_max_lon[i] = sf_max_i.lon.values
        sf_max_lat[i] = sf_max_i.lat.values

        vp_asia_i = vel_pot_mn[i,:,:]
        vp_min_i = vp_asia_i.where(vp_asia_i==vp_asia_i.min(), drop=True)
        
        vp_min_lon[i] = vp_min_i.lon.values
        vp_min_lat[i] = vp_min_i.lat.values
        
        streamfun_asia_i.plot.contourf(x='lon', y='lat', levels = np.arange(-140.,141.,5.))
        plt.plot(sf_max_lon[i], sf_max_lat[i], 'kx')
        plt.show()
        

    
    f1 = streamfun_sn[2,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(-140.,141.,10.), add_labels=False, add_colorbar=False, extend='both')
    for i in range(12):
        plt.text(sf_max_lon[i]+3, sf_max_lat[i], str(i+1))
    plt.plot(sf_max_lon, sf_max_lat, 'kx-', mew=1.5)
    plt.grid(True,linestyle=':')
    plt.colorbar(f1)
    if not land_mask==None:
        land = xr.open_dataset(land_mask)
        land.land_mask.plot.contour(x='lon', y='lat', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    
    plt.savefig(plot_dir + 'streamfun_' + run + '.pdf', format='pdf')
    plt.close()
    
    
    
    f1 = vel_pot_sn[2,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(-20.,21.,2.), add_labels=False, add_colorbar=False, extend='both')
    for i in range(12):
        plt.text(vp_min_lon[i]+3, vp_min_lat[i], str(i+1), fontsize=10)
    plt.plot(vp_min_lon, vp_min_lat, 'kx-', mew=1.5)
    plt.grid(True,linestyle=':')
    plt.colorbar(f1)
    if not land_mask==None:
        land = xr.open_dataset(land_mask)
        land.land_mask.plot.contour(x='lon', y='lat', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    
    plt.savefig(plot_dir + 'vel_pot_' + run + '.pdf', format='pdf')
    plt.close()
    
    
    

if __name__ == "__main__":
    
    #plot_sf_vp('full_qflux', land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    
    #plot_sf_vp('idealised_land', land_mask = '/scratch/rg419/Experiments/idealised_land/input/land_all.nc')
    #plot_sf_vp('idealised_land_ea_only', land_mask = '/scratch/rg419/Experiments/idealised_land/input/land_EA.nc')
    #plot_sf_vp('idealised_land_in_only', land_mask = '/scratch/rg419/Experiments/idealised_land/input/land_india.nc')
    #plot_sf_vp('idealised_land_no_af', land_mask = '/scratch/rg419/Experiments/idealised_land/input/land_asia_tibet.nc')
    #plot_sf_vp('idealised_land_no_am', land_mask = '/scratch/rg419/Experiments/idealised_land/input/land_aus_asia_tibet.nc')
    #plot_sf_vp('idealised_land_no_aus', land_mask = '/scratch/rg419/Experiments/idealised_land/input/land_africa_asia_tibet.nc')
    plot_sf_vp('idealised_land_no_sa', land_mask = '/scratch/rg419/Experiments/idealised_land/input/land_EA_tibet.nc')
                            
                            