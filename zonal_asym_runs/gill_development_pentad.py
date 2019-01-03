'''5/12/2018 Plot pentad by pentad development of the Gill pattern and with sea level pressure and precipitation
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from windspharm.xarray import VectorWind
import matplotlib.patches as patches
from data_handling_updates import model_constants as mc


def plot_gill_dev(run, land_mask=None, lev=850, qscale=100., windtype='full', ref_arrow=5, mse=False, video=False):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    data['mse'] = (mc.cp_air * data.temp + mc.L * data.sphum + mc.grav * data.height)/1000.
    data['precipitation'] = (data.precipitation*86400.)
    
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
        if mse:
            f1 = data.mse.sel(xofyear=i+1, pfull=850.).plot.contourf(x='lon', y='lat', ax=ax1, add_labels=False, add_colorbar=False, extend='both',  cmap='Blues', zorder=1, levels = np.arange(290.,341.,5.))
        else:
            f1 = data.precipitation[i,:,:].plot.contourf(x='lon', y='lat', ax=ax1, levels = np.arange(2.,21.,2.), add_labels=False, add_colorbar=False, extend='both',  cmap='Blues', zorder=1)
        #data_zanom.slp[i+4,:,:].plot.contour(x='lon', y='lat', ax=axes[i], levels = np.arange(-15.,16.,3.), add_labels=False, colors='0.5', alpha=0.5)
        ax1.contour(data_zanom.lon, data_zanom.lat, data_zanom.slp[i,:,:], levels = np.arange(0.,16.,3.), colors='0.4', alpha=0.5, zorder=2)
        ax1.contour(data_zanom.lon, data_zanom.lat, data_zanom.slp[i,:,:], levels = np.arange(-15.,0.,3.), colors='0.4', alpha=0.5, linestyle='--', zorder=2)
        if windtype=='div':
            b = ax1.quiver(data.lon[::6], data.lat[::3], uchi_zanom[i,::3,::6], vchi_zanom[i,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
            ax1.quiverkey(b, 0.,-0.5, ref_arrow, str(ref_arrow) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
        elif windtype=='rot':
            b = ax1.quiver(data.lon[::6], data.lat[::3], upsi_zanom[i,::3,::6], vpsi_zanom[i,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
            ax1.quiverkey(b, 0.,-0.5, ref_arrow, str(ref_arrow) + ' m/s', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
        elif windtype=='full':
            b = ax1.quiver(data.lon[::6], data.lat[::3], data_zanom.ucomp.sel(pfull=lev)[i,::3,::6], data_zanom.vcomp.sel(pfull=lev)[i,::3,::6], scale=qscale, angles='xy', width=0.005, headwidth=3., headlength=5., zorder=3)
            ax1.quiverkey(b, 5.,65., ref_arrow, str(ref_arrow) + ' m/s', coordinates='data', fontproperties={'weight': 'bold', 'size': 10}, color='k', labelcolor='k', labelsep=0.03, zorder=10)
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
        if mse:
            msestr = '_mse'
        if video:
            vidstr='video/'
        
        plot_dir = '/scratch/rg419/plots/zonal_asym_runs/gill_development/' + run +'/' + vidstr + windtype + msestr + '/' 
        mkdir = sh.mkdir.bake('-p')
        mkdir(plot_dir)
        
        if video:
            plt.savefig(plot_dir + 'wind_and_slp_zanom_' + str(int(data.xofyear[i])) + levstr + windtypestr + msestr + '.png', format='png')
        else:
            plt.savefig(plot_dir + 'wind_and_slp_zanom_' + str(int(data.xofyear[i])) + levstr + windtypestr + msestr + '.pdf', format='pdf')
        plt.close()


import subprocess
def make_video(filepattern, output):
    command = 'ffmpeg  -framerate 5 -y -start_number 30 -i ' + filepattern + ' -vframes 45 -c:v libx264 -r 6 -pix_fmt yuv420p -vf scale=3200:-2 ' + output    
    subprocess.call([command], shell=True)
    

if __name__ == "__main__":
    
    #plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='div', video=True)
    #plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='none', video=True)
    #plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='full', video=True)
    plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='full', mse=True, video=True)
    #plot_gill_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', windtype='full', mse=True)

    #make_video('/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/div/wind_and_slp_zanom_%02d_div.png', 
    #                  '/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/div/precip_and_slp_anom.mp4')
    
    make_video('/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/full_mse/wind_and_slp_zanom_%02d_mse.png', 
                      '/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/full_mse/mse_and_slp_anom.mp4')

    #make_video('/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/full/wind_and_slp_zanom_%02d.png', 
    #                  '/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/full/precip_and_slp_anom.mp4')
                            
                            