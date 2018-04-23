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
    
def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons

plot_dir = '/scratch/rg419/plots/flat_and_am/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

data_full = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/full_qflux.nc')
data_flat = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/flat_qflux.nc')
data_am = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/am_qflux.nc')

lons = pick_lons(data_full, [60.,150.])

data_full['totp'] = ((data_full.convection_rain + data_full.condensation_rain)*86400.).sel(lon=lons).mean('lon')
data_flat['totp'] = ((data_flat.convection_rain + data_flat.condensation_rain)*86400.).sel(lon=lons).mean('lon')
data_am['totp'] = ((data_am.convection_rain + data_am.condensation_rain)*86400.).sel(lon=lons).mean('lon')

g = 9.8
cp = 287.04/2*7
L = 2.500e6
data_full['mse'] = (data_full.temp*cp + data_full.sphum*L + data_full.height*g).mean('lon')/1000.
data_flat['mse'] = (data_flat.temp*cp + data_flat.sphum*L + data_flat.height*g).mean('lon')/1000.
data_am['mse'] = (data_am.temp*cp + data_am.sphum*L + data_am.height*g).mean('lon')/1000.

mn_dic = month_dic(1)
tickspace = range(13,72,18)
labels = [mn_dic[(k+5)/6 ] for k in tickspace]    

plevels = np.arange(2.,15.,2.)

# Begin plotting: 3 subplots

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
plt.set_cmap('RdBu_r')
#First plot
f1 = data_full.totp.plot.contourf(ax=ax1, x='xofyear', y='lat', levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
data_full.totp.plot.contour(ax=ax1, x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
cs = data_full.mse.sel(pfull=850.).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(200.,401.,20.), add_labels = False, colors='0.7', add_colorbar=False)
plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
ax1.set_ylabel('Latitude')
ax1.set_ylim(-60,60)
ax1.set_yticks(np.arange(-60.,61.,30.))
ax1.grid(True,linestyle=':')
ax1.text(-15, 60, 'a)')

f1 = data_flat.totp.plot.contourf(ax=ax2, x='xofyear', y='lat', levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
data_flat.totp.plot.contour(ax=ax2, x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
cs = data_flat.mse.sel(pfull=850.).plot.contour(ax=ax2, x='xofyear', y='lat', levels=np.arange(200.,401.,20.), add_labels = False, colors='0.7', add_colorbar=False)
plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
ax2.set_ylabel('Latitude')
ax2.set_ylim(-60,60)
ax2.set_yticks(np.arange(-60.,61.,30.))
ax2.grid(True,linestyle=':')
ax2.text(-15, 60, 'b)')

f1 = data_am.totp.plot.contourf(ax=ax3, x='xofyear', y='lat', levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
data_am.totp.plot.contour(ax=ax3, x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
cs = data_am.mse.sel(pfull=850.).plot.contour(ax=ax3, x='xofyear', y='lat', levels=np.arange(200.,401.,20.), add_labels = False, colors='0.7', add_colorbar=False)
plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
ax3.set_ylabel('Latitude')
ax3.set_ylim(-60,60)
ax3.set_yticks(np.arange(-60.,61.,30.))
ax3.grid(True,linestyle=':')
ax3.text(-15, 60, 'c)')


ax3.set_xlim((1,72))
ax3.set_xlabel('')
ax3.set_xticks(tickspace)
ax3.set_xticklabels(labels,rotation=25)
    
plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=fig.colorbar(f1, ax=(ax1, ax2, ax3), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
cb1.set_label('Precipitation, mm/day')
 
plt.savefig(plot_dir+'precip_mse_hm.pdf', format='pdf')
plt.close()        

