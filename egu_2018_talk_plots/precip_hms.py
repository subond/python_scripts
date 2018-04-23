"""
Function to plot meridional overturning, zonal wind speed, and eddy flux convergences

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from climatology import precip_mse_plot
from pylab import rcParams


rcParams['figure.figsize'] = 5, 10
rcParams['font.size'] = 14

data_ctrl = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/control_qflux.nc')
data_noam = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/no_americas.nc')
data_notp = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/no_TIP.nc')

plot_dir = '/scratch/rg419/plots/egu_2018_talk_plots/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
f1 = precip_mse_plot(data_ctrl, ax1, lonin=[60.,150.])
ax1.set_title('Control')
precip_mse_plot(data_notp, ax2, lonin=[60.,150.])
ax2.set_title('No Tibetan Plateau')
precip_mse_plot(data_noam, ax3, lonin=[60.,150.])
ax3.set_title('No America')

plt.xticks([0.,12.,24.,36.,48.,60.,72.])

plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=fig.colorbar(f1, ax=(ax1, ax2, ax3), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
cb1.set_label('Precipitation, mm/day')

plt.savefig(plot_dir+'precip_hms.pdf', format='pdf')
plt.close()     




fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
f1 = precip_mse_plot(data_ctrl, ax1, lonin=[60.,90.])
ax1.set_title('Control')
precip_mse_plot(data_notp, ax2, lonin=[60.,90.])
ax2.set_title('No Tibetan Plateau')
precip_mse_plot(data_noam, ax3, lonin=[60.,90.])
ax3.set_title('No America')

plt.xticks([0.,12.,24.,36.,48.,60.,72.])

plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=fig.colorbar(f1, ax=(ax1, ax2, ax3), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
cb1.set_label('Precipitation, mm/day')

plt.savefig(plot_dir+'precip_hms_60_90.pdf', format='pdf')
plt.close()     




fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
f1 = precip_mse_plot(data_ctrl, ax1, lonin=[90.,125.])
ax1.set_title('Control')
precip_mse_plot(data_notp, ax2, lonin=[90.,125.])
ax2.set_title('No Tibetan Plateau')
precip_mse_plot(data_noam, ax3, lonin=[90.,125.])
ax3.set_title('No America')

plt.xticks([0.,12.,24.,36.,48.,60.,72.])

plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=fig.colorbar(f1, ax=(ax1, ax2, ax3), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
cb1.set_label('Precipitation, mm/day')

plt.savefig(plot_dir+'precip_hms_90_125.pdf', format='pdf')
plt.close()     