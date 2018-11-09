"""
Load up precip and wind for each of ap_20_qflux, full_qflux, am_qflux, full and plot averages before and after onset.
1/03/2018 updated to include reference arrow.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import xarray as xr
from data_handling_updates import time_means, month_dic
import sh
from pylab import rcParams


rcParams['figure.figsize'] = 12, 7
rcParams['font.size'] = 16
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/revisions/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def precip_wind(run, ax_in, pentad):
        
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    data['rain'] = (('xofyear','lat','lon'), (data.convection_rain + data.condensation_rain)*86400.)
    
    lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= 60. and data.lon[i] < 150.]
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -30. and data.lat[i] < 60.]

    

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    
    f1 = data.rain[pentad:pentad+4,:,:].mean('xofyear').plot.contourf(x='lon', y='lat', levels = np.arange(2.,21.,2.), ax=ax_in, add_labels = False, extend='max', add_colorbar=False, cmap='Blues')
    
    land.land_mask.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_in, colors='k', add_colorbar=False, add_labels=False)
    #cs = data.t_surf[pentad:pentad+4,:,:].sel(lon=lons, lat=lats).mean('xofyear').plot.contour(x='lon', y='lat', levels = np.arange(200., 350., 10.), ax=ax_in, colors='0.6', add_labels = False, add_colorbar=False)
    #plt.clabel(cs, fontsize=12, inline_spacing=-1, fmt= '%1.0f')
    
    b = ax_in.quiver(data.sel(lon=lons, lat=lats).lon[::3], data.sel(lon=lons, lat=lats).lat[::1], data.ucomp.sel(pfull=150., lon=lons, lat=lats)[pentad:pentad+4,::1,::3].mean('xofyear'), 
                           data.vcomp.sel(pfull=150., lon=lons, lat=lats)[pentad:pentad+4,::1,::3].mean('xofyear'))#, headlength=3, headwidth=2)#, scale=500.,headwidth=5)

    ax_in.set_xlim(60,150)
    ax_in.set_ylim(-30,60)
    ax_in.set_xticks(np.arange(60.,155.,30.))
    ax_in.set_yticks(np.arange(-30.,65.,30.))
    ax_in.grid(True,linestyle=':')
    
    
    return f1, b
    
    
f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex='col', sharey='row')
    

precip_wind('ap_20_qflux', ax1, 18)
precip_wind('am_qflux', ax2, 18)
precip_wind('flat_qflux', ax3, 18)
precip_wind('full_qflux', ax4, 18)

f5, b = precip_wind('ap_20_qflux', ax5, 39)
precip_wind('am_qflux', ax6, 39)
precip_wind('flat_qflux', ax7, 39)
f1, b8 = precip_wind('full_qflux', ax8, 39)

ax5.add_patch(patches.Rectangle((60,-30), 20, 15, facecolor="white"))
ax5.quiverkey(b, 0.1,0.05,50,'50 m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03)

ax1.set_title('$ap20q$', fontsize=15)
ax2.set_title('$am20$', fontsize=15)
ax3.set_title('$flat$', fontsize=15)
ax4.set_title('$full$', fontsize=15)
    
ax1.set_ylabel('Latitude')
ax5.set_ylabel('Latitude')
ax5.set_xlabel('Longitude')
ax6.set_xlabel('Longitude')
ax7.set_xlabel('Longitude')
ax8.set_xlabel('Longitude')

ax1.text(30, 60, 'a)')
ax2.text(50, 60, 'b)')
ax3.text(50, 60, 'c)')
ax4.text(50, 60, 'd)')
ax5.text(30, 60, 'e)')
ax6.text(50, 60, 'f)')
ax7.text(50, 60, 'g)')
ax8.text(50, 60, 'h)')

    
plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.2, wspace=0.15)
#Colorbar
cb1=f.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
cb1.set_label('Precipitation, mm/day')

plot_name = plot_dir+'precip_wind_all.pdf'
plt.savefig(plot_name, format='pdf')
plt.close()