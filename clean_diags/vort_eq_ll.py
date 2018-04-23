"""
Plot absolute vorticity before and after onset at 150 hPa for full run

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means
import sh
from physics import gradients as gr
from pylab import rcParams

lev=150.
    
rcParams['figure.figsize'] = 15, 10
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/clean_diags/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

#Load in vorticity budget term means
data_vort = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_full_qflux.nc')
data_vort = data_vort * 86400.**2. #Convert to day^-2  
  
#Also load climatological data so that transient eddies can be calculated (***NB should evaluate these from daily data tomorrow***)
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/full_qflux.nc')
land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')

# Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
omega = 7.2921150e-5
f = 2 * omega * np.sin(data.lat *np.pi/180)
v_dx = gr.ddx(data.vcomp.sel(pfull=lev))  # dvdx
u_dy = gr.ddy(data.ucomp.sel(pfull=lev))  # dudy
vor = v_dx - u_dy + f
    
dvordx = gr.ddx(vor)
dvordy = gr.ddy(vor, vector=False)
    
horiz_md_mean = -86400.**2. * (data.ucomp.sel(pfull=lev) * dvordx + data.vcomp.sel(pfull=lev) * dvordy)
print 'horiz_md_mean evaluated'

div = gr.ddx(data.ucomp.sel(pfull=lev)) + gr.ddy(data.vcomp.sel(pfull=lev))
stretching_mean = -86400.**2. * vor * div
print 'stretching_mean evaluated'

transient = data_vort.horiz_md.sel(pfull=lev).values + data_vort.stretching.sel(pfull=lev).values - horiz_md_mean.values - stretching_mean.values
print 'transient evaluated'

data_vort['transient'] = (('pentad','lat','lon'), transient )	

levels = np.arange(-5.,5.1,0.5)

# Six subplots
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
plt.set_cmap('RdBu_r')
#First plot
f1 = horiz_md_mean[18:22,:,:].mean('xofyear').plot.contourf(ax=ax1, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax1.set_ylabel('Latitude')
ax1.set_ylim(-60,60)
ax1.set_xlim(60,150)
ax1.grid(True,linestyle=':')
    
#Second plot
stretching_mean[18:22,:,:].mean('xofyear').plot.contourf(ax=ax2, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax2.grid(True,linestyle=':')
ax2.set_ylim(-60,60)
ax2.set_xlim(60,150)

#Third plot
data_vort.transient[18:22,:,:].mean('pentad').plot.contourf(ax=ax3, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax3.grid(True,linestyle=':')
ax3.set_ylim(-60,60)
ax3.set_xlim(60,150)
    
#Fourth plot
horiz_md_mean[39:43,:,:].mean('xofyear').plot.contourf(ax=ax4, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax4.grid(True,linestyle=':')
ax4.set_ylabel('Latitude')
ax4.set_xlabel('Longitude')
ax4.set_ylim(-60,60)
ax4.set_xlim(60,150)
ax4.set_xticks(range(60,151,15))
    
#Fifth plot
stretching_mean[39:43,:,:].mean('xofyear').plot.contourf(ax=ax5, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax5.grid(True,linestyle=':')
ax5.set_xlabel('Longitude')
ax5.set_ylim(-60,60)
ax5.set_xlim(60,150)
ax5.set_xticks(range(60,151,15))
    
#Sixth plot
data_vort.transient[39:43,:,:].mean('pentad').plot.contourf(ax=ax6, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax6.grid(True,linestyle=':')
ax6.set_xlabel('Longitude')
ax6.set_ylim(-60,60)
ax6.set_xlim(60,150)
ax6.set_xticks(range(60,151,15))

if land:
    for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
        land.land_mask.plot.contour(x='lon', y='lat', ax=ax, colors='w', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False)
    
plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
#Colorbar
cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
cb1.set_label('day^{-2}')

figname = 'vort_budg_ll.pdf'

plt.savefig(plot_dir + figname, format='pdf')
plt.close()


