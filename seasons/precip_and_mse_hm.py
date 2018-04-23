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
    

plot_dir = '/scratch/rg419/plots/seasons/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

data_360 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')
data_720 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_2.000.nc')
data_180 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_0.500.nc')

data_360['totp'] = ((data_360.precipitation)*86400.).mean('lon')
data_720['totp'] = ((data_720.precipitation)*86400.).mean('lon')
data_180['totp'] = ((data_180.precipitation)*86400.).mean('lon')

g = 9.8
cp = 287.04/2*7
L = 2.500e6
data_360['mse'] = (data_360.temp*cp + data_360.sphum*L + data_360.height*g).mean('lon')/1000.
data_720['mse'] = (data_720.temp*cp + data_720.sphum*L + data_720.height*g).mean('lon')/1000.
data_180['mse'] = (data_180.temp*cp + data_180.sphum*L + data_180.height*g).mean('lon')/1000.

mn_dic = month_dic(1)
tickspace = range(13,72,18)
labels = [mn_dic[(k+5)/6 ] for k in tickspace]    

plevels = np.arange(2.,15.,2.)

# Begin plotting: 3 subplots

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=False)
plt.set_cmap('RdBu_r')
#First plot
f1 = data_180.totp.plot.contourf(ax=ax1, x='xofyear', y='lat', levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
data_180.totp.plot.contour(ax=ax1, x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
cs = data_180.mse.sel(pfull=850.).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(200.,401.,20.), add_labels = False, colors='0.7', add_colorbar=False)
plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
ax1.set_ylabel('Latitude')
ax1.set_ylim(-60,60)
ax1.set_yticks(np.arange(-60.,61.,30.))
ax1.grid(True,linestyle=':')
ax1.text(-15, 60, 'a)')

f1 = data_360.totp.plot.contourf(ax=ax2, x='xofyear', y='lat', levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
data_360.totp.plot.contour(ax=ax2, x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
cs = data_360.mse.sel(pfull=850.).plot.contour(ax=ax2, x='xofyear', y='lat', levels=np.arange(200.,401.,20.), add_labels = False, colors='0.7', add_colorbar=False)
plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
ax2.set_ylabel('Latitude')
ax2.set_ylim(-60,60)
ax2.set_yticks(np.arange(-60.,61.,30.))
ax2.grid(True,linestyle=':')
ax2.text(-15, 60, 'b)')

f1 = data_720.totp.plot.contourf(ax=ax3, x='xofyear', y='lat', levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
data_720.totp.plot.contour(ax=ax3, x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
cs = data_720.mse.sel(pfull=850.).plot.contour(ax=ax3, x='xofyear', y='lat', levels=np.arange(200.,401.,20.), add_labels = False, colors='0.7', add_colorbar=False)
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

