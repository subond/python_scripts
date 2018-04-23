# Load in precipitation and evaporation data and plot seasonal climatology

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh

rcParams['figure.figsize'] = 12, 10
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/paper_1_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/full_qflux.nc')
land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')

data['totp'] = (data.convection_rain + data.condensation_rain)*86400.

djf = data.sel(xofyear = range(1,13) + range(67,73)).mean('xofyear')
mam = data.sel(xofyear = range(13,31)).mean('xofyear')
jja = data.sel(xofyear = range(31,49)).mean('xofyear')
son = data.sel(xofyear = range(49,67)).mean('xofyear')

plevels = np.arange(2.,21.,2.)
# Plot precip
# Four subplots
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
#First plot
f1 = djf.totp.plot.contourf(x='lon', y='lat',ax=ax1, levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
land.land_mask.plot.contour(x='lon', y='lat',ax=ax1, levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
ax1.set_title('DJF', fontsize=17)
ax1.set_ylabel('Latitude')
ax1.set_yticks(np.arange(-90.,91.,30.))
ax1.set_xticks(np.arange(0.,361.,60.))
ax1.grid(True,linestyle=':')

#Second plot
f2 = mam.totp.plot.contourf(x='lon', y='lat',ax=ax2, levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
land.land_mask.plot.contour(x='lon', y='lat',ax=ax2, levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
ax2.set_title('MAM', fontsize=17)
ax2.set_yticks(np.arange(-90.,91.,30.))
ax2.set_xticks(np.arange(0.,361.,60.))
ax2.grid(True,linestyle=':')

#Third plot
f3 = jja.totp.plot.contourf(x='lon', y='lat',ax=ax3, levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
land.land_mask.plot.contour(x='lon', y='lat',ax=ax3, levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
ax3.set_yticks(np.arange(-90.,91.,30.))
ax3.set_xticks(np.arange(0.,361.,60.))
ax3.grid(True,linestyle=':')
ax3.set_title('JJA', fontsize=17)
ax3.set_ylabel('Latitude')
ax3.set_xlabel('Longitude')

#Fourth plot
f4 = son.totp.plot.contourf(x='lon', y='lat',ax=ax4, levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
land.land_mask.plot.contour(x='lon', y='lat',ax=ax4, levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
ax4.set_title('SON', fontsize=17)
ax4.set_yticks(np.arange(-90.,91.,30.))
ax4.set_xticks(np.arange(0.,361.,60.))
ax4.grid(True,linestyle=':')
ax4.set_xlabel('Longitude')

    
#plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=f.colorbar(f2, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.12, aspect=30, shrink=0.5)

plt.savefig(plot_dir+'precip_test.pdf', format='pdf')
plt.close()


# Plot evap
elevels = np.arange(-20.,300.,20.)
# Four subplots
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
#plt.set_cmap('RdBu_r')
#First plot
f1 = djf.flux_lhe.plot.contourf(x='lon', y='lat',ax=ax1, levels = elevels, add_colorbar=False, add_labels=False, extend='both', cmap='YlOrRd')
land.land_mask.plot.contour(x='lon', y='lat',ax=ax1, levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
ax1.set_title('DJF', fontsize=17)
ax1.set_ylabel('Latitude')
ax1.set_yticks(np.arange(-90.,91.,30.))
ax1.set_xticks(np.arange(0.,361.,60.))
ax1.grid(True,linestyle=':')

#Second plot
f2 = mam.flux_lhe.plot.contourf(x='lon', y='lat',ax=ax2, levels = elevels, add_colorbar=False, add_labels=False, extend='both', cmap='YlOrRd')
land.land_mask.plot.contour(x='lon', y='lat',ax=ax2, levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
ax2.set_title('MAM', fontsize=17)
ax2.set_yticks(np.arange(-90.,91.,30.))
ax2.set_xticks(np.arange(0.,361.,60.))
ax2.grid(True,linestyle=':')

#Third plot
f3 = jja.flux_lhe.plot.contourf(x='lon', y='lat',ax=ax3, levels = elevels, add_colorbar=False, add_labels=False, extend='both', cmap='YlOrRd')
land.land_mask.plot.contour(x='lon', y='lat',ax=ax3, levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
ax3.set_yticks(np.arange(-90.,91.,30.))
ax3.set_xticks(np.arange(0.,361.,60.))
ax3.grid(True,linestyle=':')
ax3.set_title('JJA', fontsize=17)
ax3.set_ylabel('Latitude')
ax3.set_xlabel('Longitude')

#Fourth plot
f4 = son.flux_lhe.plot.contourf(x='lon', y='lat',ax=ax4, levels = elevels, add_colorbar=False, add_labels=False, extend='both', cmap='YlOrRd')
land.land_mask.plot.contour(x='lon', y='lat',ax=ax4, levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
ax4.set_title('SON', fontsize=17)
ax4.set_yticks(np.arange(-90.,91.,30.))
ax4.set_xticks(np.arange(0.,361.,60.))
ax4.grid(True,linestyle=':')
ax4.set_xlabel('Longitude')

    
#plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=f.colorbar(f2, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.12, aspect=30, shrink=0.5)

plt.savefig(plot_dir+'evap_test.pdf', format='pdf')
plt.close()

