# 8/12/2017 Plot parts of evaporative flux to try to see where structure comes from.

from data_handling_updates import month_dic, model_constants as mc, gradients as gr
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh


g = 9.8
Rd = 287.04
Rv = 461.50
cp = Rd/2*7
L = 2.500e6
rho_cp = 1.035e3 * 3989.24495292815

def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons


def evap_parts(run, lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 6, 10
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/surface_fluxes/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    lons = pick_lons(data, lonin)
    
    # Absolute wind mag, including 1m/s vertical gust
    v_a = np.sqrt(data.ucomp_sq.sel(pfull=950.) + data.vcomp_sq.sel(pfull=950.) + 1.)
    
    # Density of lowest level air
    rho_a = 95000./Rd/data.temp.sel(pfull=950.) / (1 + (Rv - Rd)/Rd*data.sphum.sel(pfull=950.))
    
    # difference between lowest level and surface specific humidities
    q_s = Rd/Rv * 610.78/100000. * np.exp(-L/Rv * (1/data.temp.sel(pfull=950.) - 1/273.16))
    
    # Look at an average over a zonal region
    v_a = v_a.sel(lon=lons).mean('lon')
    rho_a = rho_a.sel(lon=lons).mean('lon')
    q_a = data.sphum.sel(pfull=950.).sel(lon=lons).mean('lon')
    q_s = q_s.sel(lon=lons).mean('lon')
    
    qdiff = (q_a-q_s)*1000.
    
    #(v_a*qdiff).plot.contourf(x='xofyear', y='lat', levels=np.arange(-0.05,0.05,0.005))
    #plt.show()
    #q_a.plot.contourf(x='xofyear', y='lat', levels=np.arange(0,0.02,0.001))
    #plt.figure(2)
    #q_s.plot.contourf(x='xofyear', y='lat', levels=np.arange(0,0.02,0.001))
    #plt.show()
    
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

    
    qdiff.plot.contourf(x='xofyear', y='lat', levels=np.arange(-4.5,4.6,0.5), ax=ax1, extend = 'both', add_labels=False)
    ax1.set_title('q$_{a}$ - q$_{s}$')
    
    rho_a.plot.contourf(x='xofyear', y='lat', levels=np.arange(1.,1.41,0.01), ax=ax2, extend = 'both', add_labels=False)
    ax2.set_title('rho$_{a}$')
    
    v_a.plot.contourf(x='xofyear', y='lat', levels=np.arange(0.,20.,2.),  ax=ax3, extend = 'both', add_labels=False)
    ax3.set_title('v$_{a}$')
    
    for ax in [ax1,ax2,ax3]:
        ax.set_ylim([-60,60])
        ax.set_yticks([-60,-30,0,30,60])
        ax.grid(True,linestyle=':')
    ax3.set_xticks([12,24,36,48,60,72])
    ax3.set_xlabel('Pentad')

    plt.savefig(plot_dir + 'evap_parts_' + run + '.pdf', format='pdf')
    plt.close()

    
    

if __name__ == "__main__":
    
    evap_parts('sn_1.000')