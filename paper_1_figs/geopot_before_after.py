"""
Function to plot meridional overturning, zonal wind speed, and eddy flux convergences

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from physics import gradients as gr
from pylab import rcParams

    
rcParams['figure.figsize'] = 6, 9
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
        
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/full_qflux.nc')
a= 6376.0e3 #radius used in model
coslat = np.cos(data.lat * np.pi/180)
  
height_before = 9.8 * data.height.sel(pfull=850.)[18:18+4,:,:].mean('xofyear')/1000
height_after = 9.8 * data.height.sel(pfull=850.)[39:39+4,:,:].mean('xofyear')/1000

#Geopotential gradient
dphidx_before = -86400. * 9.8 * (gr.ddx(data.height.sel(pfull=150.)))[18:18+4,:,:].mean('xofyear')
dphidx_after = -86400. * 9.8 * (gr.ddx(data.height.sel(pfull=150.)))[39:39+4,:,:].mean('xofyear')

land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')


levels = np.arange(12., 15.1, 0.25)
levels_dphidx = np.arange(-100.,100.1,20.)

# Two subplots
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
#First plot
f1 = height_before.plot.contourf(x='lon', y='lat', ax=ax1, add_labels=False, add_colorbar=False, levels=levels, extend = 'both', cmap='bwr')
land.land_mask.plot.contour(x='lon', y='lat', ax=ax1, colors='w', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False, linewidths=2)
ax1.contour(data.lon, data.lat, dphidx_before, levels = levels_dphidx, colors='k', linewidths=2)
ax1.set_ylabel('Latitude')
ax1.set_ylim(0,60)
ax1.set_xlim(60,150)
ax1.set_yticks(np.arange(0.,61.,15.))
ax1.grid(True,linestyle=':')
ax1.text(-15, 60, 'a)')

#Second plot
f2 = height_after.plot.contourf(x='lon', y='lat', ax=ax2, add_labels=False, add_colorbar=False, levels=levels, extend = 'both', cmap='bwr')
land.land_mask.plot.contour(x='lon', y='lat', ax=ax2, colors='w', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False, linewidths=2)
ax2.contour(data.lon, data.lat, dphidx_after, levels = levels_dphidx, colors='k', linewidths=2)
ax2.set_ylim(0,60)
ax2.set_xlim(60,150)
ax2.set_ylabel('Latitude')
ax2.set_xlabel('Longitude')
ax2.set_yticks(np.arange(0.,61.,15.))
ax2.set_xticks(np.arange(60.,151.,15.))
ax2.grid(True,linestyle=':')
ax2.text(-15, 60, 'b)')
    
#plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
plt.subplots_adjust(left=0.17, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=fig.colorbar(f2, ax=[ax1,ax2], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
cb1.set_label('850 hPa Geopotential, $10^3$ m$^{2}$s$^{-2}$')

plt.savefig(plot_dir+'geopotential.pdf', format='pdf')
plt.close()
    


