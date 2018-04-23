# 11/01/2018 Evaluate streamfunction, plot for JJA anomaly from time mean, and anom from zonal mean

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from physics import model_constants as mc, gradients as gr
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
    
    streamfun_sn_mean = streamfun_sn.groupby('time.month').mean('time')
    
    #time_list = [str(year) + '-08' for year in range(1979,2017)]
    #print time_list[0]
    
    #streamfun_JJA_mean = streamfun_sn.sel(time=time_list[2]).mean('time')
    #print streamfun_JJA_mean
    
    #streamfun_JJA_anom = streamfun_sn.sel(time=time_list[0]) - streamfun_JJA_mean
    #print streamfun_JJA_anom
    
    for year in range(1979,2017):

        f1 = (streamfun_sn.sel(time=str(year) + '-08') - streamfun_sn_mean.sel(month=8)).squeeze('time').plot.contourf(x='longitude', y='latitude', levels = np.arange(-20.,21.,2.), add_labels=False, add_colorbar=False, extend='both')
        plt.grid(True,linestyle=':')
        plt.colorbar(f1)

        land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    
        plt.savefig(plot_dir + 'streamfun_era_anom_' + str(year) + '.pdf', format='pdf')
        plt.close()
        
        f1 = (streamfun_sn.sel(time=str(year) + '-08') - streamfun_sn.sel(time=str(year) + '-08').mean('longitude')).squeeze('time').plot.contourf(x='longitude', y='latitude', levels = np.arange(-50.,51.,5.), add_labels=False, add_colorbar=False, extend='both')
        plt.grid(True,linestyle=':')
        plt.colorbar(f1)

        land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    
        plt.savefig(plot_dir + 'streamfun_era_zanom_' + str(year) + '.pdf', format='pdf')
        plt.close()

    

if __name__ == "__main__":
    
    plot_sf_vp(land_mask='/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc')
    