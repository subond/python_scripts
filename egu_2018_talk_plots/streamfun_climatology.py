"""
Function to plot meridional overturning, zonal wind speed, and eddy flux convergences

"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from windspharm.xarray import VectorWind

    

def plot_sf_clim(run, land_mask=None):
    
    rcParams['figure.figsize'] = 10, 8
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/egu_2018_talk_plots/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    uwnd = data.ucomp.sel(pfull=150.)
    vwnd = data.vcomp.sel(pfull=150.)
    # Create a VectorWind instance to handle the computation
    w = VectorWind(uwnd, vwnd)

    # Compute variables
    streamfun, vel_pot = w.sfvp()

    def sn_av(da):
        #Take seasonal and monthly averages
        da = da/10.**6
        da.coords['season'] = np.mod(da.xofyear + 5., 72.) // 18. 
        da_sn = da.groupby('season').mean(('xofyear'))
        return da_sn
    
    streamfun_sn = sn_av(streamfun)
    
    streamfun_sn = streamfun_sn - streamfun_sn.mean('lon')
        
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    axes = [ax1,ax2,ax3,ax4]
    title = ['DJF', 'MAM', 'JJA', 'SON']
    
    for i in range(4):
        f1 = streamfun_sn[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = np.arange(-36.,37.,3.), add_labels=False, add_colorbar=False, extend='both')
        axes[i].grid(True,linestyle=':')
        if not land_mask==None:
            land = xr.open_dataset(land_mask)
            land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
    for ax in [ax1,ax3]:
        ax.set_ylabel('Latitude')
    for ax in [ax3,ax4]:
        ax.set_xlabel('Longitude')
    
    cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    
    plt.savefig(plot_dir + 'streamfun_' + run + '.pdf', format='pdf')
    plt.close()
    

if __name__ == "__main__":
    
    plot_sf_clim('control_qflux', land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    plot_sf_clim('no_americas', land_mask = '/scratch/rg419/python_scripts/land_era/land_era_no_america.nc')
    plot_sf_clim('no_TIP', land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')