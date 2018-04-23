# Plot the surface temperature and precipitation through the year for the monsoon region

from data_handling_updates import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh

rcParams['figure.figsize'] = 8, 8
rcParams['font.size'] = 20
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')
qflux = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/ocean_qflux.nc', decode_times=False)

qflux_jja = qflux.ocean_qflux[5:8,:,:].mean('time')*(1.-land.land_mask)

# Begin plotting

f1 = qflux_jja.plot.contourf(x='lon', y='lat', levels = np.arange(-300.,301.,25.), add_colorbar=False, add_labels=False, extend='both')#, alpha=0.7)
f2 = (land.zsurf/1000).plot.contourf(x='lon', y='lat', levels = np.arange(0.,6.,0.5), add_colorbar=False, add_labels=False, extend='max', cmap='pink_r', alpha=0.8)
land.land_mask.plot.contour(x='lon', y='lat', levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
cb1=plt.colorbar(f1, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.0, aspect=30)
cb2=plt.colorbar(f2, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.2, aspect=30)

plt.ylabel('Latitude')
plt.xlabel('Longitude')
plt.xticks(np.arange(0.,361.,60.))
plt.yticks(np.arange(-90.,91.,30.))
plt.grid(True,linestyle=':')
#Colorbar
#cb1=fig.colorbar(f1, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)

plt.tight_layout()  
plt.savefig(plot_dir+'land_config.pdf', format='pdf')
plt.close()        

