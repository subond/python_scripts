# Plot up climate features (precip, winds, abs vort, tendencies) in april for ERA data

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from data_handling_updates import gradients as gr


rcParams['figure.figsize'] = 13.5, 8
rcParams['font.size'] = 18

# Load data
data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_vort_vars.nc')

land_file = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
land = xr.open_dataset( land_file)


# Calculate vorticity budget parts
stretching_mean = -1.* 86400.**2. * (data.vor * data.div)
dvordx = gr.ddx(data.vor)
dvordy = gr.ddy(data.vor, vector=False)
horiz_adv_mean = -1.* 86400.**2. * (data.ucomp * dvordx + data.vcomp * dvordy)
transient = (data.horiz_adv.values + data.stretching)* 86400.**2. - horiz_adv_mean - stretching_mean


plot_dir = '/scratch/rg419/plots/abs_vort_video/era/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


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


plt.quiver(data.lon[::3], data.lat[::1], data.ucomp[90:120,::1,::3].mean('day_of_yr'), data.vcomp[90:120,::1,::3].mean('day_of_yr'), headlength=3, headwidth=2, scale=1000., angles='xy')
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
(land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
plt.ylabel('')
plt.xlabel('')
plt.xlim(40,175)
plt.ylim(-20,60)
plt.xticks([])
plt.yticks([])
plt.tight_layout()  
plot_name = plot_dir+'wind_april.png'
plt.savefig(plot_name)
plt.close()    


data_850_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_850_clim.nc')
data_850_v = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_v/era_v_850_clim.nc')

plt.quiver(data.lon[::3], data.lat[::1], data_850_u.u_850[::1,::3,90:120].mean('day_of_yr'), data_850_v.v_850[::1,::3,90:120].mean('day_of_yr'), headlength=3, headwidth=2, scale=200., angles='xy')
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
(land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
plt.ylabel('')
plt.xlabel('')
plt.xlim(40,175)
plt.ylim(-20,60)
plt.xticks([])
plt.yticks([])
plt.tight_layout()  
plot_name = plot_dir+'wind_850_april.png'
plt.savefig(plot_name)
plt.close()



f = stretching_mean[90:120,:,:].mean('day_of_yr').plot.contourf(x='lon', y='lat', levels = np.arange(-1.5,1.75,0.25), add_labels = False, extend='both', add_colorbar=False)
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
(land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
plt.xlim(40,175)
plt.ylim(-20,60)
plt.xticks([])
plt.yticks([])
plt.tight_layout()  
plot_name = plot_dir+'stretching_april.png'
plt.savefig(plot_name)
plt.close()    


f = horiz_adv_mean[90:120,:,:].mean('day_of_yr').plot.contourf(x='lon', y='lat', levels = np.arange(-1.5,1.75,0.25), add_labels = False, extend='both', add_colorbar=False)
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
(land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
plt.xlim(40,175)
plt.ylim(-20,60)
plt.xticks([])
plt.yticks([])
plt.tight_layout()  
plot_name = plot_dir+'horiz_adv_april.png'
plt.savefig(plot_name)
plt.close()    


f = transient[90:120,:,:].mean('day_of_yr').plot.contourf(x='lon', y='lat', levels = np.arange(-1.5,1.75,0.25), add_labels = False, extend='both', add_colorbar=False)
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=5)    
(land.z[0,:,:]/9.8).plot.contourf(x='longitude', y='latitude', levels=np.arange(2500.,100001.,97500.), add_labels=False, extend='neither', add_colorbar=False, alpha=0.5, cmap='Greys_r')    
plt.xlim(40,175)
plt.ylim(-20,60)
plt.xticks([])
plt.yticks([])
plt.tight_layout()  
plot_name = plot_dir+'transient_vor_april.png'
plt.savefig(plot_name)
plt.close()




