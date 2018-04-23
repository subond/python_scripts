#Plot up horizontal advection and stretching terms in vorticity budget as a video

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from data_handling_updates import gradients as gr


rcParams['figure.figsize'] = 12, 8
rcParams['font.size'] = 16
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/vort_video/wn2/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def stretching_plot(run, ax_h, ax_s, ax_v, pentad, tibet=True):
        
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
    if tibet:
        land.zsurf.plot.contourf(x='lon', y='lat', ax=ax_h, levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')
    
    fs = stretching_mean[pentad,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(-1.5,1.6,0.5), ax=ax_s, add_labels = False, extend='both', add_colorbar=False)
    land.land_mask.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_s, colors='k', add_colorbar=False, add_labels=False)
    if tibet:
        land.zsurf.plot.contourf(x='lon', y='lat', ax=ax_s, levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')
        
    fv = (vor[pentad,:,:]*86400.).plot.contourf(x='lon', y='lat', levels = np.arange(-12.,13.,2.), ax=ax_v, add_labels = False, extend='both', add_colorbar=False)
    land.land_mask.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_v, colors='k', add_colorbar=False, add_labels=False)
    if tibet:
        land.zsurf.plot.contourf(x='lon', y='lat', ax=ax_v, levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')
    
    for ax in [ax_h, ax_s, ax_v]:
        ax.set_xlim(60,150)
        ax.set_ylim(-30,60)
        ax.set_xticks(np.arange(60.,151.,30.))
        ax.set_yticks(np.arange(-30.,61.,30.))
        ax.grid(True,linestyle=':')
        
    
    ax_h.set_xticklabels('')
    ax_s.set_xticklabels('')
    ax_v.set_xticklabels(['60','90','120','150'])
    
    
    return fh, fs, fv


for pentad in range(0,72):
    f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
    
    stretching_plot('control_qflux', ax1, ax5, ax9, pentad)
    stretching_plot('no_americas', ax2, ax6, ax10, pentad)
    #stretching_plot('frozen_am_0.4', ax3, ax7, pentad)
    [fh,fs,fv] = stretching_plot('no_TIP', ax4, ax8, ax12, pentad, tibet=False)
    
    ax1.set_title('Control', fontsize=15)
    ax2.set_title('No Am', fontsize=15)
    ax3.set_title('Frozen Am', fontsize=15)
    ax4.set_title('No TP', fontsize=15)
    
    ax1.set_ylabel('Latitude')
    ax5.set_ylabel('Latitude')
    ax5.set_xlabel('Longitude')
    ax9.set_xlabel('Longitude')
    ax10.set_xlabel('Longitude')
    ax12.set_xlabel('Longitude')
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.1, hspace=0.2, wspace=0.15)
    #Colorbar
    #cb1=f.colorbar(fh, ax=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    #cb1.set_label('Vorticity tendency, day$^{-2}$')
    
    cb1=f.colorbar(fh, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.02, aspect=20)#, shrink=0.5)
    cb1=f.colorbar(fs, ax=[ax5,ax6,ax7,ax8], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.02, aspect=20)#, shrink=0.5)
    cb1=f.colorbar(fv, ax=[ax9,ax10,ax11,ax12], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.02, aspect=20)#, shrink=0.5)
    
    plot_name = plot_dir+'vort_tends_pentad_%02d.png' % (pentad+1)
    plt.savefig(plot_name)
    plt.close()