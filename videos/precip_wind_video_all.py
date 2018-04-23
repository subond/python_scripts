import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams

rcParams['figure.figsize'] = 12, 4
rcParams['font.size'] = 16
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/precip_wind_video/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def rain_and_sst(run, ax_in, pentad, plot_land = True):
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['rain'] = (('xofyear','lat','lon'), (data.convection_rain + data.condensation_rain)*86400.)
    
    lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= 60. and data.lon[i] < 150.]
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -30. and data.lat[i] < 60.]

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
 
    f1 = data.rain[pentad,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(2.,21.,2.), ax=ax_in, add_labels = False, extend='max', add_colorbar=False, cmap='Blues')
    land.land_mask.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_in, colors='k', add_colorbar=False, add_labels=False)
    
    ax_in.quiver(data.sel(lon=lons, lat=lats).lon[::3], data.sel(lon=lons, lat=lats).lat[::1], data.ucomp.sel(pfull=150., lon=lons, lat=lats)[pentad,::1,::3], 
                           data.vcomp.sel(pfull=150., lon=lons, lat=lats)[pentad,::1,::3])#, headlength=3, headwidth=2)#,


    ax_in.set_xlim(60,150)
    ax_in.set_ylim(-30,60)
    ax_in.set_xticks(np.arange(60.,155.,30.))
    ax_in.set_yticks(np.arange(-30.,65.,30.))
    ax_in.grid(True,linestyle=':')
    
    return f1


for pentad in range(0,72):
    f, (ax1, ax2, ax3, ax4) = plt.subplots(1,4, sharey=True)
    
    rain_and_sst('ap_20_qflux', ax1, pentad)
    rain_and_sst('am_qflux', ax2, pentad)
    rain_and_sst('flat_qflux', ax3, pentad)
    f1 = rain_and_sst('full_qflux', ax4, pentad)
    
    ax1.set_title('$ap20q$', fontsize=15)
    ax2.set_title('$am20$', fontsize=15)
    ax3.set_title('$flat$', fontsize=15)
    ax4.set_title('$full$', fontsize=15)
    
    ax1.set_ylabel('Latitude')
    ax1.set_xlabel('Longitude')
    ax2.set_xlabel('Longitude')
    ax3.set_xlabel('Longitude')
    ax4.set_xlabel('Longitude')
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.9, bottom=0.1, hspace=0.2, wspace=0.15)
    #Colorbar
    cb1=f.colorbar(f1, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.2, aspect=30, shrink=0.5)
    cb1.set_label('Precipitation, mm/day')
    
    plot_name = plot_dir+'rain_wind_pentad_%02d.png' % (pentad+1)
    plt.savefig(plot_name)
    plt.close()