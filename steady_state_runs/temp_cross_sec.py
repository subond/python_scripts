"""
Plot the vertical structure of temperature at the latitude of the heating anomaly

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

ctrl_90=True

def get_temp(run, lat, ctrl_90=True):
        
    #Load data
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if ctrl_90==True:
        # Also load data from perpetual equinox run
        data_90 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_90.000.nc')
        t = data.temp.sel(lat=lat, method='nearest') - data_90.temp.sel(lat=lat, method='nearest')
        u = data.ucomp.sel(lat=lat, method='nearest') - data_90.ucomp.sel(lat=lat, method='nearest')
        w = data.omega.sel(lat=lat, method='nearest') - data_90.omega.sel(lat=lat, method='nearest')
    else:
        t = (data.temp - data.temp.mean('lon')).sel(lat=lat, method='nearest')
        u = (data.ucomp - data.ucomp.mean('lon')).sel(lat=lat, method='nearest')
        w = (data.omega - data.omega.mean('lon')).sel(lat=lat, method='nearest')
        
    return t, u, w


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

    t, u, w = get_temp(run, float(lat), ctrl=ctrl)
    
    # u is in m/s. Change to degrees longitude/day
    u = u/(mc.a * np.cos(u.lat * np.pi/180.)) * 180./np.pi * 86400.
    
    #Omega is in Pa/sec. Change to hPa/day
    w = w/100. * 86400.
    
    f1 = t.plot.contourf(x='lon', y='pfull', ax=ax, extend = 'both', add_labels=False, add_colorbar=False, levels=np.arange(-10.,11.,1.), yincrease=False)
    b = ax.quiver(u.lon[::5], w.pfull[::2], u[::2,::5], w[::2,::5], scale=1000, angles='xy')#, scale=500.,headwidth=5)
    #ax.add_patch(patches.Rectangle((0,800), 60, 200, facecolor="white"))
    #b = ax.quiverkey(b, 0.075,0.05,100,'100 hPa/day', fontproperties={'weight': 'bold', 'size': 6}, color='k', labelcolor='k', labelsep=0.03)
    #ax.set_yticks(range(-60,61,30))
    #ax.set_ylim(-60.,60.)
    ax.set_xticks(range(0,361,90))
    ax.grid(True,linestyle=':')



plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.05, hspace=0.25, wspace=0.15)

cb1=f.colorbar(f1, ax=(ax_list), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)

if ctrl_90==True:
    figname = 'temp_cross_sec.pdf'
else:
    figname = 'temp_cross_sec_zmean.pdf'
    
plt.savefig(plot_dir + figname, format='pdf')
plt.close()

