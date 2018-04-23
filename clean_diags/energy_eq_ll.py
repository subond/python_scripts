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
levels = np.arange(-0.01,0.011,0.001)

lev=850.
levels=np.arange(-0.003,0.0035,0.0005)


rcParams['figure.figsize'] = 15, 10
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/clean_diags/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

#Load in energy budget term means
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/energy_eq_full_qflux.nc')
land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')

pe_to_ke = data.pe_to_ke.sel(pfull=lev)
conv_1 = data.conv_1.sel(pfull=lev)
conv_2 = data.conv_2.sel(pfull=lev)
conv_3 = data.conv_3.sel(pfull=lev)

# Six subplots
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
plt.set_cmap('RdBu_r')

#First plot
f1 = pe_to_ke[18:22,:,:].mean('pentad').plot.contourf(ax=ax1, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax1.set_ylabel('Latitude')
ax1.set_ylim(-60,60)
ax1.set_xlim(60,150)
ax1.grid(True,linestyle=':')
    
#Second plot
conv_1[18:22,:,:].mean('pentad').plot.contourf(ax=ax2, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax2.grid(True,linestyle=':')
ax2.set_ylim(-60,60)
ax2.set_xlim(60,150)

#Third plot
conv_3[18:22,:,:].mean('pentad').plot.contourf(ax=ax3, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax3.grid(True,linestyle=':')
ax3.set_ylim(-60,60)
ax3.set_xlim(60,150)
    
#Fourth plot
pe_to_ke[39:43,:,:].mean('pentad').plot.contourf(ax=ax4, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax4.grid(True,linestyle=':')
ax4.set_ylabel('Latitude')
ax4.set_xlabel('Longitude')
ax4.set_ylim(-60,60)
ax4.set_xlim(60,150)
ax4.set_xticks(range(60,151,15))
    
#Fifth plot
conv_1[39:43,:,:].mean('pentad').plot.contourf(ax=ax5, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax5.grid(True,linestyle=':')
ax5.set_xlabel('Longitude')
ax5.set_ylim(-60,60)
ax5.set_xlim(60,150)
ax5.set_xticks(range(60,151,15))
    
#Sixth plot
conv_3[39:43,:,:].mean('pentad').plot.contourf(ax=ax6, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
ax6.grid(True,linestyle=':')
ax6.set_xlabel('Longitude')
ax6.set_ylim(-60,60)
ax6.set_xlim(60,150)
ax6.set_xticks(range(60,151,15))

for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
    land.land_mask.plot.contour(x='lon', y='lat', ax=ax, colors='w', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False)
    
plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
#Colorbar
cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
cb1.set_label('J/s')

figname = 'energy_budg_ll_' + str(int(lev)) + '.pdf'

plt.savefig(plot_dir + figname, format='pdf')
plt.close()


