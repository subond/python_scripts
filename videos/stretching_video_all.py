import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from physics import gradients as gr


rcParams['figure.figsize'] = 15, 6.5
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/stretching_video/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def stretching_plot(run, ax_in, pentad, plot_land = True):
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp.sel(pfull=150))  # dvdx
    u_dy = gr.ddy(data.ucomp.sel(pfull=150))  # dudy
    vor = v_dx - u_dy + f
    div = gr.ddx(data.ucomp.sel(pfull=150)) + gr.ddy(data.vcomp.sel(pfull=150))
    
    stretching_mean = -86400.**2. * vor * div
    
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    data['land'] = (('lat','lon'), land.land_mask)
    
    f1 = stretching_mean[pentad,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(-1.5,1.6,0.5), ax=ax_in, add_labels = False, extend='both', add_colorbar=False)
    
    if plot_land:
        data.land.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), ax=ax_in, colors='k', add_colorbar=False, add_labels=False)
    (vor*86400.)[pentad,:,:].plot.contour(x='lon', y='lat', levels = np.arange(-14.,15.,2.), ax=ax_in, colors='0.7', add_labels = False)
    
    ax_in.set_xlim(0,180)
    ax_in.set_ylim(-30,60)
    ax_in.set_xticks(np.arange(0.,185.,30.))
    ax_in.set_yticks(np.arange(-30.,65.,30.))
    
    return f1


for pentad in range(0,72):
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    
    stretching_plot('ap_2', ax1, pentad)
    stretching_plot('ap_20', ax2, pentad)
    stretching_plot('ap_20_qflux', ax3, pentad)
    stretching_plot('am_qflux', ax4, pentad)
    stretching_plot('flat_qflux', ax5, pentad)
    f1 = stretching_plot('full_qflux', ax6, pentad)
    
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
    cb1=f.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    cb1.set_label('Vorticity tendency from stretching, day$^{-2}$')
    
    plot_name = plot_dir+'stretching_pentad_%02d.png' % (pentad+1)
    plt.savefig(plot_name)
    plt.close()