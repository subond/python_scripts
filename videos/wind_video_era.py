# Make videos of era winds and of stretching, horiz adv, and transient vorticity tendencies

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from data_handling_updates import gradients as gr
from pylab import rcParams

rcParams['figure.figsize'] = 10, 5
rcParams['font.size'] = 18

# Load data
data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_vort_vars.nc')
# Take pentad means
data.coords['pentad'] = data.day_of_yr //5 + 1.  
data = data.groupby('pentad').mean(('day_of_yr')) 

land_file = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
land = xr.open_dataset( land_file)


# Calculate vorticity budget parts
stretching_mean = -1.* 86400.**2. * (data.vor * data.div)
dvordx = gr.ddx(data.vor)
dvordy = gr.ddy(data.vor, vector=False)
horiz_adv_mean = -1.* 86400.**2. * (data.ucomp * dvordx + data.vcomp * dvordy)
transient = (data.horiz_adv.values + data.stretching)* 86400.**2. - horiz_adv_mean - stretching_mean


plot_dir = '/scratch/rg419/plots/wind_video/era/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

#for i in range(0,73):
#    print i
#    plt.quiver(data.lon[::3], data.lat[::1], data.ucomp[i,::1,::3], data.vcomp[i,::1,::3], headlength=3, headwidth=2, scale=1000., angles='xy')
#    land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
#    (land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
#    plt.ylabel('Latitude')
#    plt.xlabel('Longitude')
#    plt.title('Pentad ' + str(i+1))
#    plt.xlim(0,180)
#    plt.ylim(-30,60)
#    plt.xticks(np.arange(0.,185.,30.))
#    plt.yticks(np.arange(-30.,65.,30.))
#    plt.tight_layout()  
#    plot_name = plot_dir+'wind_pentad_%02d.png' % (i+1)
#    plt.savefig(plot_name)
#    plt.close()    


plot_dir = '/scratch/rg419/plots/wind_video/era/850/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
data_850_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_850_clim.nc')
data_850_u.coords['pentad'] = data_850_u.day_of_yr //5 + 1.  
data_850_u = data_850_u.groupby('pentad').mean(('day_of_yr')) 

data_850_v = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_v/era_v_850_clim.nc')
data_850_v.coords['pentad'] = data_850_v.day_of_yr //5 + 1.  
data_850_v = data_850_v.groupby('pentad').mean(('day_of_yr')) 

for i in range(0,73):
    print i
    plt.quiver(data.lon[::3], data.lat[::1], data_850_u.u_850[i,::1,::3], data_850_v.v_850[i,::1,::3], headlength=3, headwidth=2, scale=200., angles='xy')
    land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
    (land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.title('Pentad ' + str(i+1))
    plt.xlim(0,180)
    plt.ylim(-30,60)
    plt.xticks(np.arange(0.,185.,30.))
    plt.yticks(np.arange(-30.,65.,30.))
    plt.tight_layout()  
    plot_name = plot_dir+'wind_pentad_%02d.png' % (i+1)
    plt.savefig(plot_name)
    plt.close()    


plot_dir = '/scratch/rg419/plots/stretching_video/era/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

for i in range(0,73):
    print i
    
    f = stretching_mean[i,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(-1.5,1.75,0.25), add_labels = False, extend='both', add_colorbar=False)
    land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
    (land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
    
    cb1=plt.colorbar(f)
    #cb1.set_label('Vorticity tendency from stretching, day$^{-2}$')
    
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.title('Pentad ' + str(i+1))
    plt.xlim(0,180)
    plt.ylim(-30,60)
    plt.xticks(np.arange(0.,185.,30.))
    plt.yticks(np.arange(-30.,65.,30.))
    plt.tight_layout()  
    plot_name = plot_dir+'stretching_pentad_%02d.png' % (i+1)
    plt.savefig(plot_name)
    plt.close()    


plot_dir = '/scratch/rg419/plots/horiz_adv_video/era/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

for i in range(0,73):
    print i
    f = horiz_adv_mean[i,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(-1.5,1.75,0.25), add_labels = False, extend='both', add_colorbar=False)
    land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
    (land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
    
    cb1=plt.colorbar(f)
    #cb1.set_label('Vorticity tendency from horizontal advection, day$^{-2}$')
    
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.title('Pentad ' + str(i+1))
    plt.xlim(0,180)
    plt.ylim(-30,60)
    plt.xticks(np.arange(0.,185.,30.))
    plt.yticks(np.arange(-30.,65.,30.))
    plt.tight_layout()  
    plot_name = plot_dir+'horiz_adv_pentad_%02d.png' % (i+1)
    plt.savefig(plot_name)
    plt.close()    



plot_dir = '/scratch/rg419/plots/transient_vor_video/era/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

for i in range(0,73):
    print i
    f = transient[i,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(-1.5,1.75,0.25), add_labels = False, extend='both', add_colorbar=False)
    land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
    (land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
    
    cb1=plt.colorbar(f)
    #cb1.set_label('Vorticity tendency from transient eddies, day$^{-2}$')
    
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.title('Pentad ' + str(i+1))
    plt.xlim(0,180)
    plt.ylim(-30,60)
    plt.xticks(np.arange(0.,185.,30.))
    plt.yticks(np.arange(-30.,65.,30.))
    plt.tight_layout()  
    plot_name = plot_dir+'transient_vor_pentad_%02d.png' % (i+1)
    plt.savefig(plot_name)
    plt.close()




