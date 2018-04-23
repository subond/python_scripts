"""
Plot 150 hPa wind vectors and 500 hPa vertical velocity for the anomaly experiments relative to the steady state equinox

"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr, model_constants as mc
from pylab import rcParams
import os

ctrl_90 = False

def get_wind(run, ctrl_90=True):
        
    #Load data
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if ctrl_90==True:
        # Also load data from perpetual equinox run
        data_90 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_90.000.nc')
        u_diff = data.ucomp.sel(pfull=850) - data_90.ucomp.sel(pfull=850)
        v_diff = data.vcomp.sel(pfull=850) - data_90.vcomp.sel(pfull=850)
        sst_diff  = data.t_surf - data_90.t_surf
    else:
        u_diff = data.ucomp.sel(pfull=850) - data.ucomp.sel(pfull=850).mean('lon')
        v_diff = data.vcomp.sel(pfull=850) - data.vcomp.sel(pfull=850).mean('lon')
        sst_diff  = data.t_surf - data.t_surf.mean('lon')
        
    return u_diff, v_diff, sst_diff


rcParams['figure.figsize'] = 9, 12
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/steady_state_runs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
levels = np.arange(-2.,2.1,0.25)
    
f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), 
     (ax10, ax11, ax12), (ax13, ax14, ax15), (ax16, ax17, ax18)) = plt.subplots(6, 3, sharex='col', sharey='row')
               
plt.set_cmap('RdBu_r')

ax_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17, ax18]
lat_list = [0, 5, 10, 15, 20, 25]
mag_list = [100, 200, 300]


for i in range(18):
    
    mag = mag_list[np.mod(i,3)]
    lat = lat_list[i//3]
    #print mag, lat
    ax = ax_list[i]
    
    run = 'qflux_' + str(lat) + '_' + str(mag)

    u_diff, v_diff, sst_diff = get_wind(run, ctrl_90)
    
    qflux_name = os.environ["GFDL_BASE"] + 'exp/monsoon_ss_conts/input/qflux_95_lat_' + str(lat) + '_a_' + str(mag) + '.nc'
    
    data_qflux = xr.open_dataset(qflux_name)
    
    
    f1 = sst_diff.plot.contourf(x='lon', y='lat', levels=np.arange(-10.,11.,1.), ax=ax, extend = 'both', add_labels=False, add_colorbar=False)
    data_qflux.ocean_qflux.plot.contour(x='lon', y='lat', ax=ax, add_labels=False, colors='k', levels = np.arange(50.,350.,100.))
    b = ax.quiver(u_diff.lon[::5], v_diff.lat[::2], u_diff[::2,::5], v_diff[::2,::5], scale=200., angles='xy')#,headwidth=5)
    ax.add_patch(patches.Rectangle((0,-60), 60, 20, facecolor="white"))
    b = ax.quiverkey(b, 0.075,0.05,10,'10 m/s', fontproperties={'weight': 'bold', 'size': 6}, color='k', labelcolor='k', labelsep=0.03)
    ax.set_yticks(range(-60,61,30))
    ax.set_ylim(-60.,60.)
    ax.set_xticks(range(0,361,90))
    ax.grid(True,linestyle=':')



plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.05, hspace=0.25, wspace=0.15)

cb1=f.colorbar(f1, ax=(ax_list), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)


if ctrl_90==True:
    figname = 'uv_850_sst_qflux_ll.pdf'
else:
    figname = 'uv_850_sst_qflux_ll_zmean.pdf'

plt.savefig(plot_dir + figname, format='pdf')
plt.close()

