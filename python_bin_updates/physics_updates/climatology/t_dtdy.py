# 1/12/2017 Plot hms of SST and dSSTdy with precip overlaid

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from physics import model_constants as mc, gradients as gr


g = 9.8
cp = 287.04/2*7
L = 2.500e6

def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons


def t_dtdy_plot(run, lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 6, 6
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/surface_fluxes/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    lons = pick_lons(data, lonin)
    
    try:
        precip_plot = (data.precipitation*86400.).sel(lon=lons).mean('lon')
    except:
        precip_plot = ((data.convection_rain + data.condensation_rain)*86400.).sel(lon=lons).mean('lon')
    
    data = data.sel(lon=lons).mean('lon')
    
    dTdy = gr.ddy(data.t_surf, vector=False)
    dTdy = dTdy * mc.a
    
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    
    f1 = data.t_surf.plot.contourf(x='xofyear', y='lat', levels=np.arange(270.,311.,2.5), ax=ax1, extend = 'both', add_labels=False)
    precip_plot.plot.contour(x='xofyear', y='lat', levels = np.arange(2.,15.,4.), ax=ax1, colors='k', extend = 'both', add_labels=False)
    ax1.set_title('SST')
    
    dTdy.plot.contourf(x='xofyear', y='lat', levels=np.arange(-45.,46.,5.),  ax=ax2, extend = 'both', add_labels=False)
    precip_plot.plot.contour(x='xofyear', y='lat', levels = np.arange(2.,15.,4.), ax=ax2, colors='k', extend = 'both', add_labels=False)
    ax2.set_title('dSSTdy')

    ax2.set_xlabel('Pentad')
    for ax in [ax1,ax2]:
        ax.set_ylabel('Latitude')
        ax.grid(True,linestyle=':')
    
    #plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)

    plt.savefig(plot_dir + 't_dtdy_plot_' + run + '.pdf', format='pdf')
    plt.close()


if __name__ == "__main__":
    
    t_dtdy_plot('sine_sst_10m')
    t_dtdy_plot('sn_1.000')
    

