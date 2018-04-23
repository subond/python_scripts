import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams

rcParams['figure.figsize'] = 15, 6.5
rcParams['font.size'] = 18

plot_dir = '/scratch/rg419/plots/precip_video/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def rain_and_sst(run, ax_in, pentad, plot_land = True):
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['rain'] = (('xofyear','lat','lon'), (data.convection_rain + data.condensation_rain)*86400.)

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    data['land'] = (('lat','lon'), land.land_mask)
 
    f1 = data.rain[pentad,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(2.,21.,2.), ax=ax_in, add_labels = False, extend='max', add_colorbar=False, cmap='Blues')
    if plot_land:
        data.land.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_in, colors='k', add_colorbar=False, add_labels=False)
    cs = data.t_surf[pentad,:,:].plot.contour(x='lon', y='lat', levels = np.arange(200., 350., 10.), ax=ax_in, colors='0.7', add_labels = False, add_colorbar=False)
    plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')

    ax_in.set_xlim(0,180)
    ax_in.set_ylim(-30,60)
    ax_in.set_xticks(np.arange(0.,185.,30.))
    ax_in.set_yticks(np.arange(-30.,65.,30.))
    
    return f1


for pentad in range(0,72):
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    
    rain_and_sst('ap_2', ax1, pentad)
    rain_and_sst('ap_20', ax2, pentad)
    rain_and_sst('ap_20_qflux', ax3, pentad)
    rain_and_sst('am_qflux', ax4, pentad)
    rain_and_sst('flat_qflux', ax5, pentad)
    f1 = rain_and_sst('full_qflux', ax6, pentad)
    
    ax1.set_title('ap2')
    ax2.set_title('ap20')
    ax3.set_title('ap20 + qflux')
    ax4.set_title('aquamountain')
    ax5.set_title('flat')
    ax6.set_title('full')
    
    ax1.set_ylabel('Latitude')
    ax4.set_ylabel('Latitude')
    ax4.set_xlabel('Longitude')
    ax5.set_xlabel('Longitude')
    ax6.set_xlabel('Longitude')
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.25, wspace=0.12)
    #Colorbar
    cb1=f.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    cb1.set_label('Precipitation, mm/day')
    
    plot_name = plot_dir+'rain_and_sst_pentad_%02d.png' % (pentad+1)
    plt.savefig(plot_name)
    plt.close()