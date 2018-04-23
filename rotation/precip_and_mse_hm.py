# Plot the surface temperature and precipitation through the year for the monsoon region

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh

rcParams['figure.figsize'] = 6, 10
rcParams['font.size'] = 20
rcParams['text.usetex'] = True
    

plot_dir = '/scratch/rg419/plots/rotation/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

data_10 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')
data_20 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_2.000.nc')
data_05 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_0.500.nc')

data_10['totp'] = ((data_10.precipitation)*86400.).mean('lon')
data_20['totp'] = ((data_20.precipitation)*86400.).mean('lon')
data_05['totp'] = ((data_05.precipitation)*86400.).mean('lon')

g = 9.8
cp = 287.04/2*7
L = 2.500e6
data_10['mse'] = (data_10.temp*cp + data_10.sphum*L + data_10.height*g).mean('lon')/1000.
data_20['mse'] = (data_20.temp*cp + data_20.sphum*L + data_20.height*g).mean('lon')/1000.
data_05['mse'] = (data_05.temp*cp + data_05.sphum*L + data_05.height*g).mean('lon')/1000.

mn_dic = month_dic(1)
tickspace = range(13,72,18)
labels = [mn_dic[(k+5)/6 ] for k in tickspace]    

plevels = np.arange(2.,15.,2.)

# Begin plotting: 3 subplots

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=False)
plt.set_cmap('RdBu_r')
#First plot
f1 = data_05.totp.plot.contourf(ax=ax1, x='xofyear', y='lat', levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
data_05.totp.plot.contour(ax=ax1, x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
cs = data_05.mse.sel(pfull=850.).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(200.,401.,20.), add_labels = False, colors='0.7', add_colorbar=False)
plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
ax1.set_ylabel('Latitude')
ax1.set_ylim(-60,60)
ax1.set_yticks(np.arange(-60.,61.,30.))
ax1.grid(True,linestyle=':')
ax1.text(-15, 60, 'a)')

f1 = data_10.totp.plot.contourf(ax=ax2, x='xofyear', y='lat', levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
data_10.totp.plot.contour(ax=ax2, x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
cs = data_10.mse.sel(pfull=850.).plot.contour(ax=ax2, x='xofyear', y='lat', levels=np.arange(200.,401.,20.), add_labels = False, colors='0.7', add_colorbar=False)
plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
ax2.set_ylabel('Latitude')
ax2.set_ylim(-60,60)
ax2.set_yticks(np.arange(-60.,61.,30.))
ax2.grid(True,linestyle=':')
ax2.text(-15, 60, 'b)')

f1 = data_20.totp.plot.contourf(ax=ax3, x='xofyear', y='lat', levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
data_20.totp.plot.contour(ax=ax3, x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
cs = data_20.mse.sel(pfull=850.).plot.contour(ax=ax3, x='xofyear', y='lat', levels=np.arange(200.,401.,20.), add_labels = False, colors='0.7', add_colorbar=False)
plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
ax3.set_ylabel('Latitude')
ax3.set_ylim(-60,60)
ax3.set_yticks(np.arange(-60.,61.,30.))
ax3.grid(True,linestyle=':')
ax3.text(-15, 60, 'c)')


#ax3.set_xlim((1,72))
#ax3.set_xlabel('')
#ax3.set_xticks(tickspace)
#ax3.set_xticklabels(labels,rotation=25)
    
plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=fig.colorbar(f1, ax=(ax1, ax2, ax3), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
cb1.set_label('Precipitation, mm/day')
 
plt.savefig(plot_dir+'precip_mse_hm.pdf', format='pdf')
plt.close()        

