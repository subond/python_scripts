# 11/01/2018 Evaluate streamfunction, plot for JJA with centre by month on top

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from physics import model_constants as mc, gradients as gr
from windspharm.xarray import VectorWind


month_len_list = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

def plot_sf_vp(land_mask=None):
    
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/climatology/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_mom_vars.nc')
    uwnd = data.ucomp.sel(pfull=150.)
    vwnd = data.vcomp
    # Create a VectorWind instance to handle the computation
    w = VectorWind(uwnd, vwnd)

    # Compute variables
    streamfun, vel_pot = w.sfvp()
    
    month = np.asarray([y for i in range(12) for y in month_len_list[i]*[i]])

    def mn_sn_av(da):
        #Take seasonal and monthly averages
        da = da/10.**6
        
        da.coords.update({'month': ('day_of_yr', month )})        
        da.coords['season'] = np.mod(da.month + 1., 12.) // 3. 

        da_sn = da.groupby('season').mean(('day_of_yr'))
        da_mn = da.groupby('month').mean(('day_of_yr'))
        return da_sn, da_mn
    
    streamfun_sn, streamfun_mn = mn_sn_av(streamfun)
    vel_pot_sn, vel_pot_mn = mn_sn_av(vel_pot)
        

    
    lats= [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= 5.]    
    
    sf_max_lat = np.zeros((12,))
    vp_min_lat = np.zeros((12,))
    sf_max_lon = np.zeros((12,))
    vp_min_lon = np.zeros((12,))
    
    for i in range(12):

        streamfun_asia_i = streamfun_mn[i,:,:].sel(lat=lats)
        sf_max_i = streamfun_asia_i.where(streamfun_asia_i==streamfun_asia_i.max(), drop=True)
        
        sf_max_lon[i] = sf_max_i.lon.values
        sf_max_lat[i] = sf_max_i.lat.values

        vp_asia_i = vel_pot_mn[i,:,:]
        vp_min_i = vp_asia_i.where(vp_asia_i==vp_asia_i.min(), drop=True)
        
        vp_min_lon[i] = vp_min_i.lon.values
        vp_min_lat[i] = vp_min_i.lat.values
    
    
    f1 = streamfun_sn[2,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(-140.,141.,10.), add_labels=False, add_colorbar=False, extend='both')
    for i in range(12):
        plt.text(sf_max_lon[i]+3, sf_max_lat[i], str(i+1), fontsize=10)
    plt.plot(sf_max_lon, sf_max_lat, 'kx-', mew=1.5)
    plt.grid(True,linestyle=':')
    plt.colorbar(f1)
    if not land_mask==None:
        land = xr.open_dataset(land_mask)
        land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    
    plt.savefig(plot_dir + 'streamfun_era.pdf', format='pdf')
    plt.close()
    
    
    
    f1 = vel_pot_sn[2,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(-20.,21.,2.), add_labels=False, add_colorbar=False, extend='both')
    for i in range(12):
        plt.text(vp_min_lon[i]+3, vp_min_lat[i], str(i+1), fontsize=10)
    plt.plot(vp_min_lon, vp_min_lat, 'kx-', mew=1.5)
    plt.grid(True,linestyle=':')
    plt.colorbar(f1)
    if not land_mask==None:
        land = xr.open_dataset(land_mask)
        land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    
    plt.savefig(plot_dir + 'vel_pot_era.pdf', format='pdf')
    plt.close()
    
    
    

if __name__ == "__main__":
    
    plot_sf_vp(land_mask='/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc')