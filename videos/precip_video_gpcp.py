import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams

    
    
rcParams['figure.figsize'] = 10, 5
rcParams['font.size'] = 18

plot_dir = '/scratch/rg419/plots/precip_video/gpcp/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/gpcp_detrendedpentads.nc')
land_file = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
land = xr.open_dataset( land_file)

for i in range(0,73):
    print i
    ax = data.precip_clim[i,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(2.,21.,2.), add_labels = False, extend='max', add_colorbar=False, cmap='Blues')
    land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
    (land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
    cb1=plt.colorbar(ax)
    cb1.set_label('Precipitation, mm/day')
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.title('Pentad ' + str(i+1))
    plt.xlim(0,180)
    plt.ylim(-30,60)
    plt.xticks(np.arange(0.,185.,30.))
    plt.yticks(np.arange(-30.,65.,30.))
    plt.tight_layout()  
    plot_name = plot_dir+'rain_and_sst_pentad_%02d.png' % (i+1)
    plt.savefig(plot_name)
    plt.close()    




