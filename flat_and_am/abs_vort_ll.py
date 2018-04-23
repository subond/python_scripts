"""
Plot absolute vorticity after-before onset at 150 hPa for full run

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means
from physics import mass_streamfunction
import sh
from finite_difference import cfd
from pylab import rcParams

    
rcParams['figure.figsize'] = 8, 7.3
rcParams['font.size'] = 20
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/flat_and_am/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
        
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/flat_qflux.nc')
land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')

omega = 7.2921150e-5
f = 2 * omega * np.sin(data.lat *np.pi/180)
abs_vort = (data.vor + f)*86400.

abs_vort_amb = abs_vort[39:39+4,:,:,:].mean('xofyear').sel(pfull=150) - abs_vort[18:18+4,:,:,:].mean('xofyear').sel(pfull=150)
abs_vort_after= abs_vort[39:39+4,:,:,:].mean('xofyear').sel(pfull=150)

f1 = abs_vort_amb.plot.contourf(x='lon', y='lat', extend = 'both', levels = np.arange(-5.,6.,1.), add_colorbar=False, add_labels=False)
plt.contour(data.lon, data.lat, abs_vort_after, levels = np.arange(-14.,15.,2.), colors='k', linewidths=2)
land.land_mask.plot.contour(x='lon', y='lat', levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='0.7', linewidths=2)
plt.yticks(np.arange(-90.,91.,30.))
plt.xticks(np.arange(0.,361.,60.))
plt.grid(True,linestyle=':')
plt.ylabel('Latitude')
plt.xlabel('Longitude')

#Colorbar
cb1=plt.colorbar(f1, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=40)

plt.tight_layout()

plt.savefig(plot_dir+'abs_vort_flat.pdf', format='pdf')
plt.close()
    



        
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/am_qflux.nc')
land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')

omega = 7.2921150e-5
f = 2 * omega * np.sin(data.lat *np.pi/180)
abs_vort = (data.vor + f)*86400.

abs_vort_amb = abs_vort[45:45+4,:,:,:].mean('xofyear').sel(pfull=150) - abs_vort[18:18+4,:,:,:].mean('xofyear').sel(pfull=150)
abs_vort_after= abs_vort[45:45+4,:,:,:].mean('xofyear').sel(pfull=150)

f1 = abs_vort_amb.plot.contourf(x='lon', y='lat', extend = 'both', levels = np.arange(-5.,6.,1.), add_colorbar=False, add_labels=False)
plt.contour(data.lon, data.lat, abs_vort_after, levels = np.arange(-14.,15.,2.), colors='k', linewidths=2)
land.land_mask.plot.contour(x='lon', y='lat', levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='0.7', linewidths=2)
plt.yticks(np.arange(-90.,91.,30.))
plt.xticks(np.arange(0.,361.,60.))
plt.grid(True,linestyle=':')
plt.ylabel('Latitude')
plt.xlabel('Longitude')

#Colorbar
cb1=plt.colorbar(f1, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=40)

plt.tight_layout()

plt.savefig(plot_dir+'abs_vort_am.pdf', format='pdf')
plt.close()
    




        
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ap_2.nc')
land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')

omega = 7.2921150e-5
f = 2 * omega * np.sin(data.lat *np.pi/180)
abs_vort = (data.vor + f)*86400.

abs_vort_amb = abs_vort[39:39+4,:,:,:].mean('xofyear').sel(pfull=150) - abs_vort[30:30+4,:,:,:].mean('xofyear').sel(pfull=150)
abs_vort_after= abs_vort[39:39+4,:,:,:].mean('xofyear').sel(pfull=150)

f1 = abs_vort_amb.plot.contourf(x='lon', y='lat', extend = 'both', levels = np.arange(-5.,6.,1.), add_colorbar=False, add_labels=False)
plt.contour(data.lon, data.lat, abs_vort_after, levels = np.arange(-14.,15.,2.), colors='k', linewidths=2)
land.land_mask.plot.contour(x='lon', y='lat', levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='0.7', linewidths=2)
plt.yticks(np.arange(-90.,91.,30.))
plt.xticks(np.arange(0.,361.,60.))
plt.grid(True,linestyle=':')
plt.ylabel('Latitude')
plt.xlabel('Longitude')

#Colorbar
cb1=plt.colorbar(f1, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=40)

plt.tight_layout()

plt.savefig(plot_dir+'abs_vort_ap2.pdf', format='pdf')
plt.close()
    
