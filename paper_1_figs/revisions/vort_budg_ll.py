"""
Load up vorticity budget terms for each of ap_20_qflux, full_qflux, am_qflux, full and plot averages before and after onset.

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import month_dic
import sh
from physics import gradients as gr
from pylab import rcParams


rcParams['figure.figsize'] = 12, 8
rcParams['font.size'] = 16
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/revisions/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def vort_budg(run, ax_v, ax_h, ax_s):

    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_' + run + '.nc')
    
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    
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
    
    horiz_md_hm = horiz_md_mean.sel(lon=lons).mean('lon')
    stretching_hm = stretching_mean.sel(lon=lons).mean('lon')
    abs_vort = vor.sel(lon=lons).mean('lon')*86400.
    
    abs_vort_amb = (vor[39:39+4,:,:].mean('pentad') - vor[18:18+4,:,:].mean('pentad'))*86400.
    abs_vort_after= vor[39:39+4,:,:].mean('pentad')*86400.

    levels = np.arange(-1.5,1.6,0.5)
    
    
    fv = abs_vort_amb.plot.contourf(x='lon', y='lat', extend = 'both', ax = ax_v, levels = np.arange(-5.,6.,1.), add_colorbar=False, add_labels=False)
    land.land_mask.plot.contour(x='lon', y='lat', levels = np.arange(-5.,7.,5.), ax = ax_v, add_colorbar=False, add_labels=False, colors='0.5')
    ax_v.contour(data.lon, data.lat, abs_vort_after, levels = np.arange(-20.,20.,2.), colors='k')
    ax_v.set_yticks(np.arange(-30.,61.,30.))
    ax_v.set_xticks(np.arange(60.,151.,30.))
    ax_v.set_ylim(-30,60)
    ax_v.set_xlim(60,150)
    ax_v.grid(True,linestyle=':')
    

    fh = horiz_md_mean[22:39,:,:].mean('pentad').plot.contourf(x='lon', y='lat', levels = levels, ax=ax_h, add_labels = False, extend='both', add_colorbar=False)
    land.land_mask.plot.contour(x='lon', y='lat', levels = np.arange(-5.,7.,5.), ax = ax_h, add_colorbar=False, add_labels=False, colors='0.5')
    #abs_vort.plot.contour(x='pentad', y='lat', levels = np.arange(-20.,20.,2.), ax=ax_h, add_labels = False, colors = 'k')
    
    fs = stretching_mean[22:39,:,:].mean('pentad').plot.contourf(x='lon', y='lat', levels = levels, ax=ax_s, add_labels = False, extend='both', add_colorbar=False)
    land.land_mask.plot.contour(x='lon', y='lat', levels = np.arange(-5.,7.,5.), ax = ax_s, add_colorbar=False, add_labels=False, colors='0.5')
    #abs_vort.plot.contour(x='pentad', y='lat', levels = np.arange(-20.,20.,2.), ax=ax_s, add_labels = False, colors = 'k')
  
    for ax in [ax_h, ax_s]:
        ax.set_ylim(-30,60)
        ax.set_xlim(60,150)
        ax.set_xticks(np.arange(60.,151.,30.))
        ax.set_yticks(np.arange(-30.,61.,30.))
        ax.grid(True,linestyle=':')
    
    ax_h.set_xticklabels('')
    ax_v.set_xticklabels('')
    

    return fv, fh, fs

    
f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4,  sharey='row') #sharex='col',
    

vort_budg('ap_20_qflux', ax1, ax5, ax9)
vort_budg('am_qflux', ax2, ax6, ax10)
vort_budg('flat_qflux', ax3, ax7, ax11)
[fh, fs, fv] = vort_budg('full_qflux', ax4, ax8, ax12)

ax1.set_title('$ap20q$', fontsize=15)
ax2.set_title('$am20$', fontsize=15)
ax3.set_title('$flat$', fontsize=15)
ax4.set_title('$full$', fontsize=15)

for ax in [ax1, ax5, ax9]:
    ax.set_ylabel('Latitude')

for ax in [ax9, ax10, ax11, ax12]:
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
    
plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.1, hspace=0.2, wspace=0.15)
#Colorbar
cb1=f.colorbar(fh, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.02, aspect=20)#, shrink=0.5)
cb1=f.colorbar(fs, ax=[ax5,ax6,ax7,ax8], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.02, aspect=20)#, shrink=0.5)
cb1=f.colorbar(fv, ax=[ax9,ax10,ax11,ax12], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.02, aspect=20)#, shrink=0.5)
#cb1.set_label('Vorticity tendency, day$^{-2}$')

plot_name = plot_dir+'vort_budg_fig_ll.pdf'
plt.savefig(plot_name, format='pdf')
plt.close()