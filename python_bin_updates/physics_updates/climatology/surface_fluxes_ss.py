# 1/12/2017 Plot hms of SST, precip, SW net, LW net, SH and LH

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
rho_cp = 1.035e3 * 3989.24495292815

def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons


def surface_plot(run, lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 12, 10
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/surface_fluxes/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    lons = pick_lons(data, lonin)
    
    data = data.sel(lon=lons).mean('lon')
    
    flux_lw_up = data.t_surf ** 4. * mc.stefan
    
    data.flux_sw.plot()
    
    (data.flux_lw - flux_lw_up).plot()
    
    (-1.*data.flux_t).plot()
    
    (-1.*data.flux_lhe).plot()
    
    plt.legend(['SW','LW','SENS','LH'])
    
    plt.xlabel('Latitude')
    plt.ylabel('Surface flux, W/m^2')
    plt.ylim([-300,300])
    plt.xlim([-90,90])
    
    plt.savefig(plot_dir + 'surface_fluxes_' + run + '.pdf', format='pdf')
    plt.close()



if __name__ == "__main__":
    
    surface_plot('ss_90.000')
    surface_plot('ss_91.000')
    surface_plot('ss_92.000')
    surface_plot('ss_93.000')
    surface_plot('ss_94.000')
    surface_plot('ss_95.000')
    surface_plot('ss_100.000')
    surface_plot('ss_105.000')
    surface_plot('ss_110.000')
    surface_plot('ss_115.000')



