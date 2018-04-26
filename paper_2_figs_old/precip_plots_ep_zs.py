"""
Plot precip hms for eddy permitting and zonally symmetric runs (26/03/2018)

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_mse_plot
from pylab import rcParams


data_sn1 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1_sst.nc')
data_sn1_zs = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1_sst_zs.nc')
data_sine = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sine_sst_10m.nc')
data_sine_zs = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sine_sst_10m_zs.nc')

    
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
rcParams['figure.figsize'] = 5, 10
rcParams['font.size'] = 14
    
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
f1 = precip_mse_plot(data_sn1, ax1, plot_type='precip')
ax1.set_title('sn1 SSTs')
precip_mse_plot(data_sn1_zs, ax2, plot_type='precip')
ax2.set_title('sn1 SSTs, zonally symmetric')
precip_mse_plot(data_sine, ax3, plot_type='precip')
ax3.set_title('sinusoidal SSTs')
precip_mse_plot(data_sine_zs, ax4, plot_type='precip')
ax4.set_title('sinusoidal SSTs, zonally symmetric')

for ax in [ax1, ax2, ax3, ax4]:
    ax.set_xticks([12.,24.,36.,48.,60.,72.])
    ax.set_xlim([1,72])
    ax.fill_between([32.,36.], -60, 60,
                    facecolor='k', alpha=0.3)
    ax.fill_between([39.,43.], -60, 60,
                    facecolor='k', alpha=0.3)
    ax.fill_between([46.,50.], -60, 60,
                    facecolor='k', alpha=0.3)




plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=fig.colorbar(f1, ax=(ax1, ax2, ax3, ax4), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
cb1.set_label('Precipitation, mm/day')

plt.savefig(plot_dir+'precip_mse_hm_ep_zs.pdf', format='pdf')
plt.close()        
