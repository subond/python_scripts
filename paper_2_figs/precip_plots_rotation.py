"""
Plot precip hms for rotation runs (21/02/2018)

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_mse_plot
from pylab import rcParams


data_050 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_0.500.nc')
data_075 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_0.750.nc')
data_100 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')
data_150 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_1.500.nc')
data_200 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_2.000.nc')
    
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
rcParams['figure.figsize'] = 5, 12
rcParams['font.size'] = 14
    
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True)
f1 = precip_mse_plot(data_050, ax1)
ax1.set_title('0.5xf')
precip_mse_plot(data_075, ax2)
ax2.set_title('0.75xf')
precip_mse_plot(data_100, ax3)
ax3.set_title('1xf')
precip_mse_plot(data_150, ax4)
ax4.set_title('1.5xf')
precip_mse_plot(data_200, ax5, do_xlabels=True)
ax5.set_title('2xf')
plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=fig.colorbar(f1, ax=(ax1, ax2, ax3, ax4, ax5), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
cb1.set_label('Precipitation, mm/day')

plt.savefig(plot_dir+'precip_mse_hm_rotation.pdf', format='pdf')
plt.close()        
