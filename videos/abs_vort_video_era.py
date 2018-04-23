# Make videos of era winds and of stretching, horiz adv, and transient vorticity tendencies

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams

rcParams['figure.figsize'] = 13.5, 8
rcParams['font.size'] = 18

# Load data
data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_vort_vars.nc')
# Take pentad means
data.coords['pentad'] = data.day_of_yr //5 + 1.  
data_mean = data.groupby('pentad').mean(('day_of_yr')) 

land_file = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
land = xr.open_dataset( land_file)


plot_dir = '/scratch/rg419/plots/abs_vort_video/era/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

#for i in range(0,73):
#    print i
#    f = (86400.*data_mean.vor[i,:,:]).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,2.), add_labels = False, extend='both', add_colorbar=False)
#    land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')    
#    plt.ylabel('Latitude')
#    plt.xlabel('Longitude')
#    plt.title('Pentad ' + str(i+1))
#    plt.xlim(0,180)
#    plt.ylim(-30,60)
#    plt.xticks(np.arange(0.,185.,30.))
#    plt.yticks(np.arange(-30.,65.,30.))
#    plt.tight_layout()  
#    plot_name = plot_dir+'absvor_pentad_%02d.png' % (i+1)
#    plt.savefig(plot_name)
#    plt.close()    


f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,2.), add_labels = False, extend='both', add_colorbar=False)
#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,1.), add_labels = False, extend='both', add_colorbar=False)
#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,12.1,0.1), add_labels = False, extend='both', add_colorbar=False)
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
(land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
plt.ylabel('')
plt.xlabel('')
plt.xlim(40,175)
plt.ylim(-20,60)
plt.xticks([])
plt.yticks([])
plt.tight_layout()  
plot_name = plot_dir+'absvor_may.pdf'
plt.savefig(plot_name, format='pdf')
plt.close()    


f = (86400.*data.vor[90:120,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,2.), add_labels = False, extend='both', add_colorbar=False)
#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,1.), add_labels = False, extend='both', add_colorbar=False)
#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,12.1,0.1), add_labels = False, extend='both', add_colorbar=False)
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
(land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
plt.ylabel('')
plt.xlabel('')
plt.xlim(40,175)
plt.ylim(-20,60)
plt.xticks([])
plt.yticks([])
plt.tight_layout()  
plot_name = plot_dir+'absvor_april.pdf'
plt.savefig(plot_name, format='pdf')
plt.close() 


f = (86400.*data.vor[59:90,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,2.), add_labels = False, extend='both', add_colorbar=False)
#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,1.), add_labels = False, extend='both', add_colorbar=False)
#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,12.1,0.1), add_labels = False, extend='both', add_colorbar=False)
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
(land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
plt.ylabel('')
plt.xlabel('')
plt.xlim(40,175)
plt.ylim(-20,60)
plt.xticks([])
plt.yticks([])
plt.tight_layout()  
plot_name = plot_dir+'absvor_march.pdf'
plt.savefig(plot_name, format='pdf')
plt.close()    



f = (86400.*data.vor[181:212,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,2.), add_labels = False, extend='both', add_colorbar=False)
#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,1.), add_labels = False, extend='both', add_colorbar=False)
#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,12.1,0.1), add_labels = False, extend='both', add_colorbar=False)
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
(land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
plt.ylabel('')
plt.xlabel('')
plt.xlim(40,175)
plt.ylim(-20,60)
plt.xticks([])
plt.yticks([])
plt.tight_layout()  
plot_name = plot_dir+'absvor_july.pdf'
plt.savefig(plot_name, format='pdf')
plt.close()    




#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,2.), add_labels = False, extend='both', add_colorbar=False)
#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,1.), add_labels = False, extend='both', add_colorbar=False)
#f = (86400.*data.vor[120:151,:,:].mean('day_of_yr')).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,12.1,0.1), add_labels = False, extend='both', add_colorbar=False)
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
(land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
plt.ylabel('')
plt.xlabel('')
plt.xlim(40,175)
plt.ylim(-20,60)
plt.xticks([])
plt.yticks([])
plt.tight_layout()  
plot_name = plot_dir+'asia_with_tibet_era.pdf'
plt.savefig(plot_name, format='pdf')
plt.close()    