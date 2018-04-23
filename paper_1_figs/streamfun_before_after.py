"""
Function to plot meridional overturning, zonal wind speed, and eddy flux convergences

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from physics import gradients as gr
from pylab import rcParams

    
rcParams['figure.figsize'] = 7, 9
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
        
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/full_qflux.nc')
a= 6376.0e3 #radius used in model
coslat = np.cos(data.lat * np.pi/180)
 
height_after = 9.8 * data.height[39:39+4,:,:].mean('xofyear')/1000


#sf_before = data.streamfun[18:18+4,:,:,:].mean('xofyear')
sf_after = data.streamfun[39:39+4,:,:,:].mean('xofyear')

land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')

levels_low = np.arange(-2.5e7, 2.6e7, 0.5e7)
levels_high = np.arange(-1.e8, 1.e8, 0.2e8)
levels = np.arange(12., 15.1, 0.25)
levels_phi_high = np.arange(120.,150.1,2.)

# Two subplots
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
#First plot
f1 = height_after.sel(pfull=850.).plot.contourf(x='lon', y='lat', ax=ax1, add_labels=False, levels=levels, extend = 'both', cmap='bwr')
land.land_mask.plot.contour(x='lon', y='lat', ax=ax1, colors='0.7', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False, linewidths=2)
ax1.contour(data.lon, data.lat, height_after.sel(pfull=150.), levels = levels_phi_high, colors='k', linewidths=2)
#f1 = sf_before.sel(pfull=850.).plot.contourf(x='lon', y='lat', ax=ax1, add_labels=False, add_colorbar=False, extend = 'both', cmap='bwr', levels=levels_low)
#land.land_mask.plot.contour(x='lon', y='lat', ax=ax1, colors='w', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False, linewidths=2)
#ax1.contour(data.lon, data.lat, sf_before.sel(pfull=150.), colors='k', linewidths=2, levels=levels_high)
ax1.set_ylabel('Latitude')
ax1.set_xticks(np.arange(0.,361.,60.))
ax1.set_yticks(np.arange(-90.,91.,30.))
ax1.grid(True,linestyle=':')
ax1.text(-50, 90, 'a)')
ax1.set_title('Geopotential', fontsize=17)

#Second plot
f2 = sf_after.sel(pfull=850.).plot.contourf(x='lon', y='lat', ax=ax2, add_labels=False, extend = 'both', cmap='bwr', levels=levels_low)
land.land_mask.plot.contour(x='lon', y='lat', ax=ax2, colors='0.7', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False, linewidths=2)
ax2.contour(data.lon, data.lat, sf_after.sel(pfull=150.), colors='k', linewidths=2, levels=levels_high)
ax2.set_ylabel('Latitude')
ax2.set_xlabel('Longitude')
ax2.set_xticks(np.arange(0.,361.,60.))
ax2.set_yticks(np.arange(-90.,91.,30.))
ax2.grid(True,linestyle=':')
ax2.text(-50, 90, 'b)')
ax2.set_title('Horizontal streamfunction', fontsize=17)

    
#plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
plt.subplots_adjust(left=0.17, right=0.95, top=0.95, bottom=0.1, hspace=0.2)
#Colorbar
#cb1=fig.colorbar(f2, ax=[ax1,ax2], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
#cb1.set_label('850 hPa Geopotential, $10^3$ m$^{2}$s$^{-2}$')

plt.savefig(plot_dir+'horiz_streamfun.pdf', format='pdf')
plt.close()
    
