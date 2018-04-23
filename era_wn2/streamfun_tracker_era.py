# 11/01/2018 Evaluate streamfunction, plot for JJA with centre by month on top

#from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
#from physics import model_constants as mc, gradients as gr
from windspharm.xarray import VectorWind
import pandas as pd



def plot_sf_vp(land_mask=None):
    
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/era_wn2/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    land = xr.open_dataset(land_mask)
    
    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_200_mm.nc')
    print data.time.units
    uwnd = data.u.load().squeeze('level')
    uwnd = uwnd[0:456,:,:]
    data_sn = data.resample(time='Q-NOV').mean()
    uwnd_sn = data_sn.u.load().squeeze('level')
    
    
    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_v/era_v_200_mm.nc')
    vwnd = data.v.load()
    data_sn = data.resample(time='Q-NOV').mean()
    vwnd_sn = data_sn.v.load()
    

    # Create a VectorWind instance to handle the computation
    w = VectorWind(uwnd, vwnd)
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    streamfun = streamfun/10.**6
    
    # Create a VectorWind instance to handle the computation
    w = VectorWind(uwnd_sn, vwnd_sn)
    # Compute variables
    streamfun_sn, vel_pot_sn = w.sfvp()
    streamfun_sn = streamfun_sn/10.**6
    
    lats = [data.latitude[i] for i in range(len(data.latitude)) if data.latitude[i] >= 5.]    
    
    
    for year in range(1979,2017):
        
        sf_max_lat = np.zeros((12,))
        sf_max_lon = np.zeros((12,))
        
        for month in range(1,13):
            streamfun_asia_i = streamfun.sel(time=str(year)+'-%02d' % month , latitude=lats)
            
            sf_max_i = streamfun_asia_i.where(streamfun_asia_i==streamfun_asia_i.max(), drop=True)
        
            sf_max_lon[month-1] = sf_max_i.longitude.values
            sf_max_lat[month-1] = sf_max_i.latitude.values
        
        f1 = streamfun_sn.sel(time=str(year) + '-08').squeeze('time').plot.contourf(x='longitude', y='latitude', levels = np.arange(-140.,141.,10.), add_labels=False, add_colorbar=False, extend='both')
        for i in range(12):
            plt.text(sf_max_lon[i]+3, sf_max_lat[i], str(i+1), fontsize=10)
        plt.plot(sf_max_lon, sf_max_lat, 'kx-', mew=1.5)
        plt.grid(True,linestyle=':')
        plt.colorbar(f1)

        land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    
        plt.savefig(plot_dir + 'streamfun_era_' + str(year) + '.pdf', format='pdf')
        plt.close()

    

if __name__ == "__main__":
    
    plot_sf_vp(land_mask='/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc')
    