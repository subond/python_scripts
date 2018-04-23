"""
Plot absolute vorticity before and after onset at 150 hPa for full run

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means
from physics import mass_streamfunction
import sh
from finite_difference import cfd
from pylab import rcParams

    
rcParams['figure.figsize'] = 6, 9
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/clean_diags/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
        
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/full_qflux.nc')
land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')

omega = 7.2921150e-5
f = 2 * omega * np.sin(data.lat *np.pi/180)
abs_vort = (data.vor + f)*86400.

abs_vort_before = abs_vort[18:18+4,:,:,:].mean('xofyear').sel(pfull=150)
abs_vort_after= abs_vort[39:39+4,:,:,:].mean('xofyear').sel(pfull=150)


# Two subplots
f, (ax1, ax2) = plt.subplots(2, sharex=True)
plt.set_cmap('RdBu_r')
#First plot
f1 = ax1.contourf(data.lon, data.lat, abs_vort_before, extend = 'both', levels = np.arange(-14.,15.,2.))
land.land_mask.plot.contour(ax=ax1, x='lon', y='lat', levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
ax1.set_ylabel('Latitude')
ax1.set_ylim(-60,60)
ax1.grid(True,linestyle=':')

#Third plot
f2 = ax2.contourf(data.lon, data.lat, abs_vort_after, extend = 'both', levels = np.arange(-14.,15.,2.))
land.land_mask.plot.contour(ax=ax2, x='lon', y='lat', levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
ax2.grid(True,linestyle=':')
ax2.invert_yaxis()
#ax2.set_xlim(0,150)
ax2.set_ylim(-60,60)
ax2.set_ylabel('Latitude')
ax2.set_xlabel('Longitude')

    
#plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=f.colorbar(f2, ax=[ax1,ax2], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
cb1.set_label('Absolute vorticity')

plt.savefig(plot_dir+'abs_vort_full.pdf', format='pdf')
plt.close()
    


