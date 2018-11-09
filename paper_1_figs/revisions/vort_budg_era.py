"""
Load variables for vorticity budget and save in one dataset for ease of use
"""
import numpy as np
import xarray as xr
import sh
from data_handling_updates import gradients as gr
from pylab import rcParams
import matplotlib.pyplot as plt

    
rcParams['figure.figsize'] = 15, 8
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/paper_1_figs/revisions/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
filename = '/scratch/rg419/obs_and_reanalysis/era_vort_vars.nc'
data = xr.open_dataset(filename, decode_times=False)

# Take pentad means
data.coords['pentad'] = data.day_of_yr //5 + 1.  
data = data.groupby('pentad').mean(('day_of_yr')) 

# Choose lon range
#lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= 60. and data.lon[i] < 150.]
#lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= 30. and data.lon[i] < 60.]
lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= 15. and data.lon[i] < 60.]
    
# Calculate time mean stretching
stretching_mean = -1.* 86400.**2. * (data.vor * data.div)
stretching_hm = stretching_mean.sel(lon=lons).mean('lon')

dvordx = gr.ddx(data.vor)
dvordy = gr.ddy(data.vor, vector=False)

horiz_adv_mean = -1.* 86400.**2. * (data.ucomp * dvordx + data.vcomp * dvordy)
horiz_adv_hm = horiz_adv_mean.sel(lon=lons).mean('lon')

transient_s_hm = (86400.**2. * data.stretching.values - stretching_mean).sel(lon=lons).mean('lon')
transient_h_hm = (86400.**2. * data.horiz_adv.values - horiz_adv_mean).sel(lon=lons).mean('lon')

transient_hm = transient_s_hm + transient_h_hm    

    
tickspace = [13, 31, 49, 68]
labels = ['Mar', 'Jun', 'Sep', 'Dec']    
    
levels = np.arange(-1.5,1.6,0.25)
# Four subplots
f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex='col', sharey='row')
plt.set_cmap('RdBu_r')
#First plot
f2=horiz_adv_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax1, extend = 'both', add_labels=False)
ax1.set_ylabel('Latitude')
ax1.set_ylim(-60,60)
ax1.set_title('Horizontal advection', fontsize=17)
ax1.set_yticks(np.arange(-60,61,30))
ax1.grid(True,linestyle=':')
ax1.text(-15, 60, 'a)')

#Second plot
dvordy.sel(lon=lons).mean('lon').plot.contourf(x='pentad', y='lat', levels=np.arange(-1.e-10,1.05e-10,0.1e-10), ax=ax2, extend = 'both', add_labels=False)
ax2.set_yticks(np.arange(-60,61,30))
ax2.set_title('$\partial (\overline{\zeta} + f) /\partial y$', fontsize=17)
ax2.grid(True,linestyle=':')
ax2.text(-7, 60, 'b)')
    

#Third plot
data.vcomp.sel(lon=lons).mean('lon').plot.contourf(x='pentad', y='lat', levels=np.arange(-12.,13,2.), ax=ax3, extend = 'both', add_labels=False)
ax3.set_ylim(-60,60)
ax3.set_yticks(np.arange(-60,61,30))
ax3.set_title('$\overline{v}$', fontsize=17)
ax3.grid(True,linestyle=':')
ax3.text(-7, 60, 'c)')
    
#Fourth plot
stretching_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax4, extend = 'both', add_labels=False)
ax4.set_ylabel('Latitude')
ax4.set_ylim(-60,60)
ax4.set_yticks(np.arange(-60,61,30))
ax4.set_xticks(tickspace)
ax4.set_xticklabels(labels,rotation=25)
ax4.set_title('Vortex stretching', fontsize=17)
ax4.grid(True,linestyle=':')
ax4.text(-15, 60, 'd)')
    
#Fifth plot
(data.div.sel(lon=lons).mean('lon')*86400.).plot.contourf(x='pentad', y='lat', levels=np.arange(-0.6,0.7,0.1), ax=ax5, extend = 'both', add_labels=False)
ax5.set_ylim(-60,60)
ax5.set_yticks(np.arange(-60,61,30))
ax5.set_xticks(tickspace)
ax5.set_xticklabels(labels,rotation=25)
ax5.grid(True,linestyle=':')
ax5.set_title('Divergence', fontsize=17)
ax5.text(-7, 60, 'e)')
    
#Sixth plot
(data.vor.sel(lon=lons).mean('lon')*86400.).plot.contourf(x='pentad', y='lat', levels=np.arange(-12.,13.,2.), ax=ax6, extend = 'both', add_labels=False)
ax6.set_ylim(-60,60)
ax6.set_yticks(np.arange(-60,61,30))
ax6.set_xticks(tickspace)
ax6.set_xticklabels(labels,rotation=25)
ax6.grid(True,linestyle=':')
ax6.set_title('Absolute vorticity', fontsize=17)
ax6.text(-7, 60, 'f)')

    
#Seventh plot
transient_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax7, extend = 'both', add_labels=False)
#ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
ax7.set_ylabel('Latitude')
ax7.set_ylim(-60,60)
ax7.set_yticks(np.arange(-60,61,30))
ax7.set_xticks(tickspace)
ax7.set_xticklabels(labels,rotation=25)
ax7.set_title('Transient eddy vorticity tendency', fontsize=17)
ax7.grid(True,linestyle=':')
ax7.text(-15, 60, 'g)')

#Seventh plot
transient_h_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax8, extend = 'both', add_labels=False)
#ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
ax8.set_ylim(-60,60)
ax8.set_yticks(np.arange(-60,61,30))
ax8.set_xticks(tickspace)
ax8.set_xticklabels(labels,rotation=25)
ax8.set_title('Transient horizontal adv.', fontsize=17)
ax8.grid(True,linestyle=':')
ax8.text(-7, 60, 'h)')

#Seventh plot
transient_s_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax9, extend = 'both', add_labels=False)
#ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
ax9.set_ylim(-60,60)
ax9.set_yticks(np.arange(-60,61,30))
ax9.set_xticks(tickspace)
ax9.set_xticklabels(labels,rotation=25)
ax9.set_title('Transient vortex stretching', fontsize=17)
ax9.grid(True,linestyle=':')
ax9.text(-7, 60, 'i)')


plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.1, hspace=0.25, wspace=0.15)
    
figname = 'vort_breakdown_era_southafrica.pdf'

plt.savefig(plot_dir + figname, format='pdf')
plt.close()    