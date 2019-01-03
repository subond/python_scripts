'''6/12/2018 Plot pentad by pentad development of half_shallow Rossby no and precip
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from windspharm.xarray import VectorWind
import matplotlib.patches as patches
from data_handling_updates import gradients as gr, model_constants as mc


def plot_vort_dev(run, land_mask=None, lev=200, video=False, threed=True):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['precipitation'] = (data.precipitation*86400.)
    
    zon_adv = data.ucomp * gr.ddx(data.ucomp)
    merid_adv = data.vcomp * gr.ddy(data.ucomp)
    vert_adv = data.omega * gr.ddp(data.ucomp)
    
    sinphi = np.sin(data.lat * np.pi/180.)
    f = 2.* mc.omega*sinphi
    if threed:
        rossby =  (zon_adv + merid_adv + vert_adv)/(f*data.vcomp)
    else:
        rossby = merid_adv/(f*data.vcomp)
    # Start figure with 1 subplots
    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 14

    for i in range(72):
        fig, ax1 = plt.subplots()
        title = 'Pentad ' + str(int(data.xofyear[i]))
        
        f1 = rossby.sel(xofyear=i+1, pfull=lev).plot.contourf(x='lon', y='lat', ax=ax1, add_labels=False, add_colorbar=False, extend='both', zorder=1, levels = np.arange(0.,1.1,0.1))
        
        data.precipitation.sel(xofyear=i+1).plot.contour(x='lon', y='lat', ax=ax1, add_labels=False, extend='both', zorder=1, levels=np.arange(2.,21.,2.), colors='k')        
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
    
        vidstr=''
        if video:
            vidstr='video/'
        
        if threed:
            plot_dir = '/scratch/rg419/plots/zonal_asym_runs/gill_development/' + run +'/' + vidstr + '/rossby3d/' 
        else:
            plot_dir = '/scratch/rg419/plots/zonal_asym_runs/gill_development/' + run +'/' + vidstr + '/rossby/' 
        mkdir = sh.mkdir.bake('-p')
        mkdir(plot_dir)
        
        if video:
            plt.savefig(plot_dir + 'rossby_and_precip_' + str(int(data.xofyear[i])) + '.png', format='png')
        else:
            plt.savefig(plot_dir + 'rossby_and_precip_' + str(int(data.xofyear[i])) + '.pdf', format='pdf')
        plt.close()


import subprocess
def make_video(filepattern, output):
    command = 'ffmpeg  -framerate 5 -y -start_number 30 -i ' + filepattern + ' -vframes 45 -c:v libx264 -r 6 -pix_fmt yuv420p -vf scale=3200:-2 ' + output    
    subprocess.call([command], shell=True)
    

if __name__ == "__main__":

    plot_vort_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', threed=True)
    plot_vort_dev('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', video=True, threed=True)

    make_video('/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/rossby3d/rossby_and_precip_%02d.png', 
                      '/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/rossby3d/rossby_and_precip.mp4')

    
                         