''' 
6/11/2017
Plot the peak overturning at 500 hPa as a function of time
'''

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import mass_streamfunction, get_edge_psi, overturning_plot


# control
data_10 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')

# mlds
data_2 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ap_2.nc')
data_20 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ap_20.nc')

# seasons
data_sn05 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_0.500.nc')
data_sn20 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_2.000.nc')

# rotation
data_rt05 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_0.500.nc')
data_rt20 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_2.000.nc')

# newer runs
data_zs = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/zs_sst.nc')
data_co2 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ap10_co2.nc')
data_qflux = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ap10_qflux.nc')


plot_dir = '/scratch/rg419/plots/itcz_plots/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    

rcParams['font.size'] = 18
rcParams['text.usetex'] = True


def cell_strength_plot(psi_max, plotname, clrs=None, key_labels=None):
    n = len(psi_max)
    if clrs==None:
        clrs = ['k'] * n
    
    for i in range(n):
        psi_max[i].plot(color=clrs[i])
    
    if not key_labels==None:
        plt.legend(key_labels, loc='upper left')
    
    plt.ylabel('Peak mass streamfunction, 10$^9$ kg/s')
    plt.xlabel('Pentad')

    
    plt.savefig(plot_dir+'max_overturning_' + plotname + '.pdf', format='pdf')
    plt.close()


def get_max_overturning(data):
    n = len(data)
    
    psi_max=[]
    for i in range(n):
        vars = get_edge_psi(data[i], lev=500.)
        psi_max.append(vars[1])
    return psi_max

psi_max = get_max_overturning([data_2,data_10,data_20])
cell_strength_plot(psi_max, 'mlds', clrs=['b','k','r'], key_labels=['2m', '10m', '20m'])

psi_max = get_max_overturning([data_sn05,data_10,data_sn20])
cell_strength_plot(psi_max, 'sns', clrs=['b','k','r'], key_labels=['180 day', '360 day', '720 day'])

psi_max = get_max_overturning([data_rt05,data_10,data_rt20])
cell_strength_plot(psi_max, 'rts', clrs=['b','k','r'], key_labels=['0.5*$\Omega$', '1*$\Omega$', '2*$\Omega$'])

psi_max = get_max_overturning([data_10, data_zs, data_co2, data_qflux])
cell_strength_plot(psi_max, 'other', clrs=['b','k','r', 'g'], key_labels=['10m', 'axisymmetric', '4*CO2', 'Q-flux'])

