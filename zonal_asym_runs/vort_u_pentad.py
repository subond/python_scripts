'''6/12/2018 Plot pentad by pentad development of half_shallow upper level vorticity and wind vectors
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from windspharm.xarray import VectorWind
import matplotlib.patches as patches
from data_handling_updates import gradients as gr, model_constants as mc


def plot_vort_dev(run, land_mask=None, lev=200, qscale=150., windtype='full', ref_arrow=10., video=False):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    sinphi = np.sin(data.lat * np.pi/180.)
    zeta = (2.* mc.omega*sinphi -1.* gr.ddy(data.ucomp)) * 86400.
    
    # Take zonal anomaly
    data_zanom = data - data.mean('lon')
    
    # Get rotational and divergent components of the flow
    w = VectorWind(data.ucomp.sel(pfull=lev), data.vcomp.sel(pfull=lev))
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    uchi_zanom = (uchi - uchi.mean('lon')).sortby('lat')
    vchi_zanom = (vchi - vchi.mean('lon')).sortby('lat')
    upsi_zanom = (upsi - upsi.mean('lon')).sortby('lat')
    vpsi_zanom = (vpsi - vpsi.mean('lon')).sortby('lat')
    
    # Start figure with 1 subplots
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 14

    for i in range(72):
        fig, ax1 = plt.subplots()
        title = 'Pentad ' + str(int(data.xofyear[i]))
        
        f1 = zeta.sel(xofyear=i+1, pfull=lev).plot.contourf(x='lon', y='lat', ax=ax1, add_labels=False, add_colorbar=False, extend='both', zorder=1, levels = np.arange(-10.,10.,2.))
        
        if windtype=='div':
            b = ax1.quiver(data.lon[::6], data.lat[::3], uchi_zanom[i,::3,::6], vchi_zanom[i,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
            ax1.quiverkey(b, 0.,-0.5, ref_arrow, str(ref_arrow) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
        elif windtype=='rot':
            b = ax1.quiver(data.lon[::6], data.lat[::3], upsi_zanom[i,::3,::6], vpsi_zanom[i,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
            ax1.quiverkey(b, 0.,-0.5, ref_arrow, str(ref_arrow) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
        elif windtype=='full':
            b = ax1.quiver(data.lon[::6], data.lat[::3], data_zanom.ucomp.sel(pfull=lev)[i,::3,::6], data_zanom.vcomp.sel(pfull=lev)[i,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
            ax1.quiverkey(b, 0.,-0.5, ref_arrow, str(ref_arrow) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
        else:
            windtype='none'
        ax1.grid(True,linestyle=':')
        ax1.set_ylim(-60.,60.)
        ax1.set_yticks(np.arange(-60.,61.,30.))
        ax1.set_xticks(np.arange(0.,361.,90.))
        ax1.set_title(title)
        if not land_mask==None:
            land = xr.open_dataset(land_mask)
            land.land_mask.plot.contour(x='lon', y='lat', ax=ax1, levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')    
        ax1.set_ylabel('Latitude')
        ax1.set_xlabel('Longitude')
    
        plt.subplots_adjust(left=0.1, right=0.97, top=0.93, bottom=0.05, hspace=0.25, wspace=0.2)
        cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.15, aspect=60, shrink=0.5)
    
        levstr=''; windtypestr=''; msestr=''; vidstr=''
        if lev != 850:
            levstr = '_' + str(lev)
        if windtype != 'full':
            windtypestr = '_' + windtype
        if video:
            vidstr='video/'
        
        plot_dir = '/scratch/rg419/plots/zonal_asym_runs/gill_development/' + run +'/' + vidstr + windtype + '/' 
        mkdir = sh.mkdir.bake('-p')
        mkdir(plot_dir)
        
        if video:
            plt.savefig(plot_dir + 'wind_and_vort_zanom_' + str(int(data.xofyear[i])) + levstr + windtypestr + '.png', format='png')
        else:
            plt.savefig(plot_dir + 'wind_and_vort_zanom_' + str(int(data.xofyear[i])) + levstr + windtypestr + '.pdf', format='pdf')
        plt.close()


import subprocess
def make_video(filepattern, output):
    command = 'ffmpeg  -framerate 5 -y -start_number 30 -i ' + filepattern + ' -vframes 45 -c:v libx264 -r 6 -pix_fmt yuv420p -vf scale=3200:-2 ' + output    
    subprocess.call([command], shell=True)
    

if __name__ == "__main__":

    #plot_vort_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='full')
    #plot_vort_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='full', video=True)
    #plot_vort_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='div')
    #plot_vort_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='div', video=True)
    #plot_vort_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='rot')
    #plot_vort_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='rot', video=True)

    #make_video('/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/full/wind_and_vort_zanom_%02d_200.png', 
    #                  '/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/full/wind_and_vort_200.mp4')

    make_video('/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/div/wind_and_vort_zanom_%02d_200_div.png', 
                      '/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/div/wind_and_vort_200.mp4')
    
    make_video('/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/rot/wind_and_vort_zanom_%02d_200_rot.png', 
                      '/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/rot/wind_and_vort_200.mp4')
    
                         