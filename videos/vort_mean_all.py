#Plot up horizontal advection and stretching terms in vorticity budget as a video

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from physics import gradients as gr


rcParams['figure.figsize'] = 12, 8
rcParams['font.size'] = 16
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/vort_video/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def stretching_plot(run, ax_h, ax_s, pentad):
        
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
    
    levels = np.arange(-1.5,1.6,0.5)
    
    fh = horiz_md_mean[pentad,:,:].plot.contourf(x='lon', y='lat', levels = levels, ax=ax_h, add_labels = False, extend='both', add_colorbar=False)
    land.land_mask.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_h, colors='k', add_colorbar=False, add_labels=False)
    
    fs = stretching_mean[pentad,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(-1.5,1.6,0.5), ax=ax_s, add_labels = False, extend='both', add_colorbar=False)
    land.land_mask.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_s, colors='k', add_colorbar=False, add_labels=False)
    
    for ax in [ax_h, ax_s]:
        ax.set_xlim(60,150)
        ax.set_ylim(-30,60)
        ax.set_xticks(np.arange(60.,151.,30.))
        ax.set_yticks(np.arange(-30.,61.,30.))
        ax.grid(True,linestyle=':')
        
    
    ax_h.set_xticklabels('')
    ax_s.set_xticklabels(['60','90','120','150'])
    
    
    return fh, fs


for pentad in range(0,72):
    f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex='col', sharey='row')
    
    stretching_plot('ap_20_qflux', ax1, ax5, pentad)
    stretching_plot('am_qflux', ax2, ax6, pentad)
    stretching_plot('flat_qflux', ax3, ax7, pentad)
    [fh,fs] = stretching_plot('full_qflux', ax4, ax8, pentad)
    
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
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.1, hspace=0.2, wspace=0.15)
    #Colorbar
    cb1=f.colorbar(fh, ax=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    cb1.set_label('Vorticity tendency, day$^{-2}$')
    
    plot_name = plot_dir+'stretching_pentad_%02d.png' % (pentad+1)
    plt.savefig(plot_name)
    plt.close()