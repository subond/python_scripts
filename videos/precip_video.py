import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams

    
def rain_and_sst(run, plot_land = True):
    
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 18
    
    plot_dir = '/scratch/rg419/plots/precip_video/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['rain'] = (('xofyear','lat','lon'), (data.convection_rain + data.condensation_rain)*86400.)

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    data['land'] = (('lat','lon'), land.land_mask)
 
    for i in range(0,72):
        print i
        ax = data.rain[i,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(2.,15.,2.), add_labels = False, extend='max', add_colorbar=False, cmap='Blues')
        if plot_land:
            data.land.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), colors='k', add_colorbar=False, add_labels=False)
        cs = data.t_surf[i,:,:].plot.contour(x='lon', y='lat', levels = np.arange(200., 350., 10.), colors='0.7', add_labels = False, add_colorbar=False)
        plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
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
        plot_name = plot_dir+'rain_and_sst_pentad_%02d.png' % i
        plt.savefig(plot_name)
        plt.close()    


#rain_and_sst('full_qflux')
rain_and_sst('ap_20_qflux')
rain_and_sst('flat_qflux')
rain_and_sst('am_qflux')
rain_and_sst('ap_2')
rain_and_sst('ap_20')
