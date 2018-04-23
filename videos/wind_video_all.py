import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams

rcParams['figure.figsize'] = 15, 6.5
rcParams['font.size'] = 18

plot_dir = '/scratch/rg419/plots/wind_video/850/'
plot_dir = '/scratch/rg419/plots/wind_video/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def wind(run, ax_in, pentad, plot_land = True):
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    data['land'] = (('lat','lon'), land.land_mask)
    
    ax_in.quiver(data.lon[::3], data.lat[::1], data.ucomp.sel(pfull=150.)[pentad,::1,::3], data.vcomp.sel(pfull=150.)[pentad,::1,::3], headlength=3, headwidth=2, scale=500., angles='xy')
    if plot_land:
        data.land.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_in, colors='k', add_colorbar=False, add_labels=False)

    ax_in.set_xlim(0,180)
    ax_in.set_ylim(-30,60)
    ax_in.set_xticks(np.arange(0.,185.,30.))
    ax_in.set_yticks(np.arange(-30.,65.,30.))
    
    #return f1


for pentad in range(0,72):
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    
    wind('ap_2', ax1, pentad)
    wind('ap_20', ax2, pentad)
    wind('ap_20_qflux', ax3, pentad)
    wind('am_qflux', ax4, pentad)
    wind('flat_qflux', ax5, pentad)
    wind('full_qflux', ax6, pentad)
    
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
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.1, hspace=0.25, wspace=0.12)
    #Colorbar
    #cb1=f.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    #cb1.set_label('Precipitation, mm/day')
    
    plot_name = plot_dir+'wind_pentad_%02d.png' % (pentad+1)
    plt.savefig(plot_name)
    plt.close()