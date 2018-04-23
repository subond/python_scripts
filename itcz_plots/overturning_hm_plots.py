''' 
6/11/2017
Plot the overturning at 500 hPa, along with the peak overturning location, and the cell edge diagnostic
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
data_sine = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sine_sst_10m.nc')
data_sine_zs = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sine_sst_10m_zs.nc')

data_dry_ep = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/dry_ep.nc')
data_dry_zs = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/dry_zs.nc')


plot_dir = '/scratch/rg419/plots/itcz_plots/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    

rcParams['font.size'] = 18
rcParams['text.usetex'] = True


def plot_multiple(data, plotname, lonin=[-1.,361.], do_xlabels=False, month_labels=True, sharex=True, thresh=0.1, nh=False):
    n = len(data)
    rcParams['figure.figsize'] = 6, 10/3*n
    
    fig, axes = plt.subplots(n, sharex=sharex)
    axes = axes.tolist()
    
    f1 = overturning_plot(data[0], axes[0], lonin=lonin, do_xlabels=do_xlabels, month_labels=month_labels, thresh=thresh, nh=nh)
    for i in range(1, n-1):
        overturning_plot(data[i], axes[i], lonin=lonin, do_xlabels=do_xlabels, month_labels=month_labels, thresh=thresh, nh=nh)
    overturning_plot(data[n-1], axes[n-1], lonin=lonin, do_xlabels=True, month_labels=month_labels, thresh=thresh, nh=nh)
    
    plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
    #Colorbar
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
    
    cb1.set_label('Mass streamfunction, 10$^9$ kg/s')

    
    plt.savefig(plot_dir+'overturning_hm_' + plotname + '.pdf', format='pdf')
    plt.close()
    


#data = [data_2, data_10, data_20]
#plot_multiple(data, 'mlds')

#data = [data_sn05, data_10, data_sn20]
#plot_multiple(data, 'sns', sharex=False, month_labels=False)

#data = [data_rt05, data_10, data_rt20]
#plot_multiple(data, 'rts')

#data = [data_10, data_zs, data_co2, data_qflux]
#plot_multiple(data, 'other')

data = [data_10, data_zs, data_sine, data_sine_zs]
plot_multiple(data, 'symmetry')

#data = [data_dry_ep, data_dry_zs]
#plot_multiple(data, 'dry')
