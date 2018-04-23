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
        u_diff = data.ucomp.sel(pfull=150) - data_90.ucomp.sel(pfull=150)
        v_diff = data.vcomp.sel(pfull=150) - data_90.vcomp.sel(pfull=150)
        omega_diff  = (data.omega.sel(pfull=500) - data_90.omega.sel(pfull=500))*86400./100. # Convert to hPa/day
    else:
        u_diff = data.ucomp.sel(pfull=150) - data.ucomp.sel(pfull=150).mean('lon')
        v_diff = data.vcomp.sel(pfull=150) - data.vcomp.sel(pfull=150).mean('lon')
        omega_diff  = (data.omega.sel(pfull=500) - data.omega.sel(pfull=150).mean('lon'))*86400./100. # Convert to hPa/day
        
    return u_diff, v_diff, omega_diff



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

    u_diff, v_diff, omega_diff = get_wind(run, ctrl_90=ctrl_90)
    
    qflux_name = os.environ["GFDL_BASE"] + 'exp/monsoon_ss_conts/input/qflux_95_lat_' + str(lat) + '_a_' + str(mag) + '.nc'
    
    data_qflux = xr.open_dataset(qflux_name)
    
    
    f1 = omega_diff.plot.contourf(x='lon', y='lat', levels=np.arange(-120.,121.,10.), ax=ax, extend = 'both', add_labels=False, add_colorbar=False)
    data_qflux.ocean_qflux.plot.contour(x='lon', y='lat', ax=ax, add_labels=False, colors='k', levels = np.arange(50.,350.,100.))
    if ctrl_90=='True':
        b = ax.quiver(u_diff.lon[::5], v_diff.lat[::2], u_diff[::2,::5], v_diff[::2,::5], scale=1000., angles='xy')#, scale=500.,headwidth=5)
        ax.add_patch(patches.Rectangle((0,-60), 60, 20, facecolor="white"))
        b = ax.quiverkey(b, 0.075,0.05,50,'50 m/s', fontproperties={'weight': 'bold', 'size': 6}, color='k', labelcolor='k', labelsep=0.03)
    else:
        b = ax.quiver(u_diff.lon[::5], v_diff.lat[::2], u_diff[::2,::5], v_diff[::2,::5], scale=500., angles='xy')#, scale=500.,headwidth=5)
        ax.add_patch(patches.Rectangle((0,-60), 60, 20, facecolor="white"))
        b = ax.quiverkey(b, 0.075,0.05,25,'25 m/s', fontproperties={'weight': 'bold', 'size': 6}, color='k', labelcolor='k', labelsep=0.03)
    ax.set_yticks(range(-60,61,30))
    ax.set_ylim(-60.,60.)
    ax.set_xticks(range(0,361,90))
    ax.grid(True,linestyle=':')



plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.05, hspace=0.25, wspace=0.15)

cb1=f.colorbar(f1, ax=(ax_list), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)

if ctrl_90==True:
    figname = 'uv_150_w_500_qflux_ll.pdf'
else:
    figname = 'uv_150_w_500_qflux_ll_zmean.pdf'


plt.savefig(plot_dir + figname, format='pdf')
plt.close()

