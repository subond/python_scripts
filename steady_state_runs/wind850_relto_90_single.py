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

def get_wind(run):
        
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


rcParams['figure.figsize'] = 9, 6
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/steady_state_runs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
levels = np.arange(-2.,2.1,0.25)
    
plt.set_cmap('RdBu_r')

def plot_wind(run, qflux_tag):
    
    u_diff, v_diff, sst_diff = get_wind(run)
    
    qflux_name = os.environ["GFDL_BASE"] + 'exp/monsoon_ss_conts/input/qflux_95_lat_' + qflux_tag + '.nc'
    
    data_qflux = xr.open_dataset(qflux_name)
    
    f1 = sst_diff.plot.contourf(x='lon', y='lat', levels=np.arange(-1.,1.1,0.1), extend = 'both', add_labels=False, add_colorbar=False)
    data_qflux.ocean_qflux.plot.contour(x='lon', y='lat', add_labels=False, colors='k', levels = np.arange(50.,350.,100.))
    b = plt.quiver(u_diff.lon[::2], v_diff.lat[::1], u_diff[::1,::2], v_diff[::1,::2], scale=10., angles='xy')#,headwidth=5)
    #plt.add_patch(patches.Rectangle((0,-60), 60, 20, facecolor="white"))
    b = plt.quiverkey(b, 0.075,0.05,1,'1 m/s', fontproperties={'weight': 'bold', 'size': 6}, color='k', labelcolor='k', labelsep=0.03)
    plt.yticks(range(-60,61,30))
    plt.ylim(-60.,60.)
    plt.xticks(range(0,361,90))
    plt.grid(True,linestyle=':')
    
    cb1=plt.colorbar(f1, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)
    
    if ctrl_90==True:
        figname = 'uv_850_sst_' + run + '.pdf'
    else:
        figname = 'uv_850_sst_' + run + '_zmean.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()


plot_wind('qflux_0_200_small', '0_a_200_small')