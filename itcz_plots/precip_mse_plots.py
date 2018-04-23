''' 
6/11/2017
Plot precipitation and/or MSE
'''

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from climatology import precip_mse_plot


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


plot_dir = '/scratch/rg419/plots/itcz_plots/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    

rcParams['font.size'] = 18

var_dic = {'lonin': [-1.,361.], 
          'plot_type': None, 
          'precip_contour': 8., 
          'p_cent': True, 
          'mse_max': True, 
          'month_labels': True}


def plot_multiple(data, plotname, do_xlabels=False, sharex=True, var_dic = var_dic):
    n = len(data)
    rcParams['figure.figsize'] = 6, 10/3*n
    
    fig, axes = plt.subplots(n, sharex=sharex)
    axes = axes.tolist()

    f1 = precip_mse_plot(data[0], axes[0], **var_dic)
    for i in range(1, n-1):
        precip_mse_plot(data[i], axes[i], **var_dic)
    precip_mse_plot(data[n-1], axes[n-1], do_xlabels=True, **var_dic)
    
    plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
    #Colorbar
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
    
    if var_dic['plot_type']=='mse':
        cb1.set_label('Moist static energy, kJ/kg')
    else:
        cb1.set_label('Precipitation, mm/day')
    
    if var_dic['plot_type']==None:
        plt.savefig(plot_dir+'precip_mse_hm_' + plotname + '.pdf', format='pdf')
    elif var_dic['plot_type']=='precip':
        plt.savefig(plot_dir+'precip_hm_' + plotname + '.pdf', format='pdf')
    elif var_dic['plot_type']=='mse':
        plt.savefig(plot_dir+'mse_hm_' + plotname + '.pdf', format='pdf')
     
    plt.close()
    
for plot_type in [None, 'precip',  'mse']:
    var_dic['plot_type'] = plot_type
    
    #data = [data_2, data_10, data_20]
    #plot_multiple(data, 'mlds')
    
    #var_dic['month_labels'] = False
    #data = [data_sn05, data_10, data_sn20]
    #plot_multiple(data, 'sns', sharex=False)
    
    #var_dic['month_labels'] = True
    #data = [data_rt05, data_10, data_rt20]
    #plot_multiple(data, 'rts')
    
    #data = [data_10, data_zs, data_co2, data_qflux]
    #plot_multiple(data, 'other')

    data = [data_10, data_zs, data_sine, data_sine_zs]
    plot_multiple(data, 'symmetry')



