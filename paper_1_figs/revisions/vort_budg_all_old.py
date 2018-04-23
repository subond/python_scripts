"""
Load up vorticity budget terms for each of ap_20_qflux, full_qflux, am_qflux, full and plot averages before and after onset.

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams


rcParams['figure.figsize'] = 12, 14
rcParams['font.size'] = 16
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/revisions/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def vort_budg(run, ax_s, ax_h, pentad):
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['rain'] = (('xofyear','lat','lon'), (data.convection_rain + data.condensation_rain)*86400.)
    
    lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= 60. and data.lon[i] < 150.]
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -30. and data.lat[i] < 60.]

    
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp.sel(pfull=150))  # dvdx
    u_dy = gr.ddy(data.ucomp.sel(pfull=150))  # dudy
    vor = v_dx - u_dy + f
    div = gr.ddx(data.ucomp.sel(pfull=150)) + gr.ddy(data.vcomp.sel(pfull=150))
    
    dvordx = gr.ddx(vor)
    dvordy = gr.ddy(vor, vector=False)
    
    horiz_md_mean = -86400.**2. * (data.ucomp.sel(pfull=150) * dvordx + data.vcomp.sel(pfull=150) * dvordy)
    stretching_mean = -86400.**2. * vor * div
    

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    
    f1 = stretching_mean[pentad:pentad+4,:,:].mean('xofyear').plot.contourf(x='lon', y='lat', levels = np.arange(-5.,5.1,1.), ax=ax_s, add_labels = False, extend='both', add_colorbar=False)
    
    land.land_mask.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_s, colors='k', add_colorbar=False, add_labels=False)

    cs = (vor*86400.)[pentad:pentad+4,:,:].sel(lon=lons, lat=lats).mean('xofyear').plot.contour(x='lon', y='lat', levels = np.arange(-14.,15.,2.), ax=ax_s, colors='0.6', add_labels = False, add_colorbar=False)
    plt.clabel(cs, fontsize=12, inline_spacing=-1, fmt= '%1.0f')
    
    
    f1 = horiz_md_mean[pentad:pentad+4,:,:].mean('xofyear').plot.contourf(x='lon', y='lat', levels = np.arange(-5.,5.1,1.), ax=ax_h, add_labels = False, extend='both', add_colorbar=False)
    
    land.land_mask.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_h, colors='k', add_colorbar=False, add_labels=False)

    cs = (vor*86400.)[pentad:pentad+4,:,:].sel(lon=lons, lat=lats).mean('xofyear').plot.contour(x='lon', y='lat', levels = np.arange(-14.,15.,2.), ax=ax_h, colors='0.6', add_labels = False, add_colorbar=False)
    plt.clabel(cs, fontsize=12, inline_spacing=-1, fmt= '%1.0f')
    
    for ax in [ax_s, ax_h]:
        ax.set_xlim(60,150)
        ax.set_ylim(-30,60)
        ax.set_xticks(np.arange(60.,155.,30.))
        ax.set_yticks(np.arange(-30.,65.,30.))
    
    return f1
    
    
f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16)) = plt.subplots(4, 4, sharex='col', sharey='row')
    

vort_budg('ap_20_qflux', ax1, ax9, 18)
vort_budg('am_qflux', ax2, ax10, 18)
vort_budg('flat_qflux', ax3, ax11, 18)
vort_budg('full_qflux', ax4, ax12, 18)

vort_budg('ap_20_qflux', ax5, ax13, 39)
vort_budg('am_qflux', ax6, ax14, 39)
vort_budg('flat_qflux', ax7, ax15, 39)
f1 = vort_budg('full_qflux', ax8, ax16, 39)

ax1.set_title('$ap20 + qflux$', fontsize=15)
ax2.set_title('$aquamountain$', fontsize=15)
ax3.set_title('$flat$', fontsize=15)
ax4.set_title('$full$', fontsize=15)

for ax in [ax1, ax5, ax9, ax13]:
    ax.set_ylabel('Latitude')

for ax in [ax13,ax14,ax15,ax16]:
    ax.set_xlabel('Longitude')


ax1.text(30, 60, 'a)')
ax2.text(50, 60, 'b)')
ax3.text(50, 60, 'c)')
ax4.text(50, 60, 'd)')
ax5.text(30, 60, 'e)')
ax6.text(50, 60, 'f)')
ax7.text(50, 60, 'g)')
ax8.text(50, 60, 'h)')
ax9.text(30, 60, 'i)')
ax10.text(50, 60, 'j)')
ax11.text(50, 60, 'k)')
ax12.text(50, 60, 'l)')
ax13.text(30, 60, 'm)')
ax14.text(50, 60, 'n)')
ax15.text(50, 60, 'o)')
ax16.text(50, 60, 'p)')
    
plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.25, wspace=0.15)
#Colorbar
cb1=f.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8, ax9,ax10,ax11,ax12,ax13,ax14,ax15,ax16], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.05, aspect=30, shrink=0.5)
cb1.set_label('Vorticity tendency, day$^{-2}$')

plot_name = plot_dir+'vort_budg_all.pdf'
plt.savefig(plot_name, format='pdf')
plt.close()