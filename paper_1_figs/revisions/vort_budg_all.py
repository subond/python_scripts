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


def vort_budg(run, ax_h, ax_s, ax_v):

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

    levels = np.arange(-1.5,1.6,0.25)
    
    
    fh = horiz_md_hm.plot.contourf(x='pentad', y='lat', levels = levels, ax=ax_h, add_labels = False, extend='both', add_colorbar=False)
    abs_vort.plot.contour(x='pentad', y='lat', levels = np.arange(-20.,20.,2.), ax=ax_h, add_labels = False, colors = 'k')
    
    fs = stretching_hm.plot.contourf(x='pentad', y='lat', levels = levels, ax=ax_s, add_labels = False, extend='both', add_colorbar=False)
    abs_vort.plot.contour(x='pentad', y='lat', levels = np.arange(-20.,20.,2.), ax=ax_s, add_labels = False, colors = 'k')
  
    for ax in [ax_h, ax_s]:
        ax.set_ylim(-60,60)
        ax.set_xticks(range(13,72,18))
        ax.set_yticks(np.arange(-60.,61.,30.))
        ax.grid(True,linestyle=':')
    
    ax_h.set_xticklabels('')
    
    fv = abs_vort_amb.plot.contourf(x='lon', y='lat', extend = 'both', ax = ax_v, levels = np.arange(-5.,6.,1.), add_colorbar=False, add_labels=False)
    land.land_mask.plot.contour(x='lon', y='lat', levels = np.arange(-5.,7.,5.), ax = ax_v, add_colorbar=False, add_labels=False, colors='0.5')
    ax_v.contour(data.lon, data.lat, abs_vort_after, levels = np.arange(-20.,20.,2.), colors='k')
    ax_v.set_yticks(np.arange(-90.,91.,30.))
    ax_v.set_xticks(np.arange(0.,181.,30.))
    ax_v.set_ylim(-60,60)
    ax_v.set_xlim(0,180)
    ax_v.grid(True,linestyle=':')
        
        
    return fh, fs, fv

    
f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3,  sharey='row') #sharex='col',
    

vort_budg('ap_20_qflux', ax1, ax4, ax7)
vort_budg('am_qflux', ax2, ax5, ax8)
[fh, fs, fv] = vort_budg('flat_qflux', ax3, ax6, ax9)

ax1.set_title('$ap20 + qflux$', fontsize=15)
ax2.set_title('$aquamountain$', fontsize=15)
ax3.set_title('$flat$', fontsize=15)

for ax in [ax1, ax4, ax7]:
    ax.set_ylabel('Latitude')

for ax in [ax7, ax8, ax9]:
    ax.set_xlabel('Longitude')

for ax in [ax4, ax5, ax6]:
    mn_dic = month_dic(1)
    labels = [mn_dic[(k+5)/6 ] for k in range(13,72,18)]
    ax.set_xticklabels(labels,rotation=25)

ax1.text(-15, 60, 'a)')
ax2.text(-5, 60, 'b)')
ax3.text(-5, 60, 'c)')
ax4.text(-15, 60, 'd)')
ax5.text(-5, 60, 'e)')
ax6.text(-5, 60, 'f)')
ax7.text(-40, 60, 'g)')
ax8.text(-15, 60, 'h)')
ax9.text(-15, 60, 'i)')

    
plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.1, hspace=0.25, wspace=0.15)
#Colorbar
cb1=f.colorbar(fh, ax=[ax1,ax2,ax3], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.02, aspect=20)#, shrink=0.5)
cb1=f.colorbar(fs, ax=[ax4,ax5,ax6], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.02, aspect=20)#, shrink=0.5)
cb1=f.colorbar(fv, ax=[ax7,ax8,ax9], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.02, aspect=20)#, shrink=0.5)
#cb1.set_label('Vorticity tendency, day$^{-2}$')

plot_name = plot_dir+'vort_budg_fig.pdf'
plt.savefig(plot_name, format='pdf')
plt.close()