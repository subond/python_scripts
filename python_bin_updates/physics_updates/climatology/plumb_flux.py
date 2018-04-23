# 7/02/2018 Evaluate the Plumb flux for a given simulation

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from data_handling_updates import model_constants as mc, gradients as gr


def plumb_flux(data, plev=200.):
    
    prefac_a = plev/1000. * np.cos(np.pi*data.lat/180.)
    prefac_b = 1/(2*mc.omega*mc.a*np.sin(np.pi*data.lat/180.))
    
    u_anom = (data.ucomp - data.ucomp.mean('lon')).sel(pfull=plev)
    v_anom = (data.vcomp - data.vcomp.mean('lon')).sel(pfull=plev)
    z_anom = mc.grav*(data.height - data.height.mean('lon')).sel(pfull=plev)
    
    f_x = prefac_a * (v_anom**2 - prefac_b * gr.ddx(v_anom * z_anom))
    
    f_y = prefac_a * (-u_anom * v_anom + prefac_b * gr.ddx(u_anom * z_anom))
    
    return f_x, f_y


def plot_plumb_flux(run, ax_in, plev=200., time_av = 'pentad', time_index=41, land_mask=None):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    if time_av == 'season':
        # Take seasonal averages
        data.coords['season'] = np.mod(data.xofyear + 5., 72.) // 18. 
        #print data.season
        data = data.groupby('season').mean(('xofyear'))
        if time_index > 3:
            time_index = 2   # Plot summer as default 
    
    f_x, f_y = plumb_flux(data, plev=plev)
    
    ax_in.quiver(data.lon[::5], data.lat[::2], f_x[::2,time_index,::5], f_y[::2,time_index,::5], angles='xy', scale=100.)
    if not land_mask==None:
        land = xr.open_dataset(land_mask)
        land.land_mask.plot.contour(x='lon', y='lat', ax=ax_in, levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
    
    

def plot_plumb_flux_diff(run1, run2, ax_in, plev=200., time_av = 'pentad', time_index=41, land_mask=None):
    
    data1 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run1 + '.nc')
    data2 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run2 + '.nc')
    
    if time_av == 'season':
        # Take seasonal averages
        data1.coords['season'] = np.mod(data1.xofyear + 5., 72.) // 18. 
        data2.coords['season'] = np.mod(data2.xofyear + 5., 72.) // 18. 
        #print data.season
        data1 = data1.groupby('season').mean(('xofyear'))
        data2 = data2.groupby('season').mean(('xofyear'))
        if time_index > 3:
            time_index = 2   # Plot summer as default
            
    f_x1, f_y1 = plumb_flux(data1, plev=plev)
    f_x2, f_y2 = plumb_flux(data2, plev=plev)
    
    f_x = f_x2 - f_x1
    f_y = f_y2 - f_y1
    
    ax_in.quiver(data1.lon, data1.lat, f_x[:,time_index,:], f_y[:,time_index,:], angles='xy', scale=100.)
    if not land_mask==None:
        land = xr.open_dataset(land_mask)
        land.land_mask.plot.contour(x='lon', y='lat', ax=ax_in, levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        


def plumb_seasons_plot(run1, run2=None, plev=200., land_mask=None):
    
    rcParams['figure.figsize'] = 10, 8
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/climatology/plumb_flux/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    axes = [ax1,ax2,ax3,ax4]
    title = ['DJF', 'MAM', 'JJA', 'SON']
    
    if run2 == None:
        for i in range(4):
            plot_plumb_flux(run1, ax_in=axes[i], plev=plev, time_av='season', time_index=i, land_mask=land_mask)
            axes[i].set_title(title[i])
        plt.savefig(plot_dir + 'plumb_flux_' + run1 + '.pdf', format='pdf')
        plt.close()
    else:
        for i in range(4):
            plot_plumb_flux_diff(run1, run2, ax_in=axes[i], plev=plev, time_av='season', time_index=i, land_mask=land_mask)
            axes[i].set_title(title[i])
        plt.savefig(plot_dir + 'plumb_flux_' + run1 + '_vs_' + run2 + '.pdf', format='pdf')
        plt.close()
    
    

if __name__ == "__main__":
    
    #plot_plumb_flux('idealised_1hill', plev=850., time_av='season', land_mask = '/scratch/rg419/Experiments/wavenumber_2/input/1hill.nc')
    #plot_plumb_flux('idealised_2hill', plev=850., time_av='season', land_mask = '/scratch/rg419/Experiments/wavenumber_2/input/2hill.nc')
    #plot_plumb_flux_diff('idealised_1hill', 'idealised_2hill', plev=850., time_av='season', land_mask = '/scratch/rg419/Experiments/wavenumber_2/input/2hill.nc')
    
    plumb_seasons_plot('idealised_1hill', land_mask = '/scratch/rg419/Experiments/wavenumber_2/input/1hill.nc', plev=850)
    plumb_seasons_plot('idealised_2hill', land_mask = '/scratch/rg419/Experiments/wavenumber_2/input/2hill.nc', plev=850)
    plumb_seasons_plot('idealised_1hill', run2='idealised_2hill', land_mask = '/scratch/rg419/Experiments/wavenumber_2/input/2hill.nc', plev=850)