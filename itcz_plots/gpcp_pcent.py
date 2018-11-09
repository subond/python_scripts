''' 
12/07/2018 Plot the precipitation centroid in lat-lon on GPCP data at selected pentads
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from climatology import precip_centroid_ll
from data_handling_updates import cell_area_from_xar, gradients as gr

data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/gpcp_detrendedpentads.nc')
data['precipitation'] = data.precip_clim
data['latb'] = np.arange(-90.,90.5)
data['lonb'] = np.arange(-180.,180.5)

area = cell_area_from_xar(data, radius = 6371.0e3)[0]
data['area'] = xr.DataArray(area, coords=[data.lat.values, data.lon.values], dims=['lat', 'lon'])



plot_dir = '/scratch/rg419/plots/itcz_plots/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

land = xr.open_dataset('/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
land_2 = np.zeros((64,128))
land_2[:,0:64] = land.land_mask[:,64:]
land_2[:,64:] = land.land_mask[:,0:64]
lon_2 = np.zeros((128,))
lon_2[0:64] = land.lon[64:]-360.
lon_2[64:] = land.lon[0:64]

land_mask = xr.DataArray(land_2, coords=[land.lat.values, lon_2], dims=['lat', 'lon'])

data = precip_centroid_ll(data, lat_bound=30.)

dpcentdt_india = gr.ddt(data.p_cent.sel(lon=np.arange(70.5,90.5)).mean('lon')) * 86400. 
dpcentdt_sasia = gr.ddt(data.p_cent.sel(lon=np.arange(90.5,130.5)).mean('lon')) * 86400. 
#dpcentdt_sasia = gr.ddt(data.p_cent.sel(lon=np.arange(110.5,120.5)).mean('lon')) * 86400. 
dpcentdt_sa = gr.ddt(data.p_cent.sel(lon=np.arange(-70.5,-45.5)).mean('lon')) * 86400. 
#dpcentdt_africa = gr.ddt(data.p_cent.sel(lon=np.arange(17.5,35.5)).mean('lon')) * 86400. 
dpcentdt_africa = gr.ddt(data.p_cent.sel(lon=np.arange(15.5,60.5)).mean('lon')) * 86400. 


data_2p5 = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/mld_2.5.nc')
data_10 = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/sn_1.000.nc')
data_15 = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/mld_15.nc')
data_20 = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/mld_20.nc')

data_2p5 = precip_centroid_ll(data_2p5, lat_bound=30.)
data_10 = precip_centroid_ll(data_10, lat_bound=30.)
data_15 = precip_centroid_ll(data_15, lat_bound=30.)
data_20 = precip_centroid_ll(data_20, lat_bound=30.)

dpcentdt_2p5 = gr.ddt(data_2p5.p_cent.mean('lon')) * 86400. 
dpcentdt_10 = gr.ddt(data_10.p_cent.mean('lon')) * 86400. 
dpcentdt_15 = gr.ddt(data_15.p_cent.mean('lon')) * 86400. 
dpcentdt_20 = gr.ddt(data_20.p_cent.mean('lon')) * 86400. 


# Set figure parameters
rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 14
# Start figure with 4 subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    
ax1.plot(data.p_cent.sel(lon=np.arange(70.5,90.5)).mean('lon'), dpcentdt_india, 'kx', mew=2, ms=10, alpha=0.5)
ax1.plot(data.p_cent.sel(lon=np.arange(70.5,90.5)).mean('lon'), dpcentdt_india, alpha=0.3, color='k', linewidth=2)
ax1.plot(data_2p5.p_cent.mean('lon'), dpcentdt_2p5, linewidth=2)
ax1.plot(data_10.p_cent.mean('lon'), dpcentdt_10, linewidth=2)
ax1.plot(data_15.p_cent.mean('lon'), dpcentdt_15, linewidth=2)
ax1.plot(data_20.p_cent.mean('lon'), dpcentdt_20, linewidth=2)
ax1.set_title('India (70-90E)')

ax2.plot(data.p_cent.sel(lon=np.arange(90.5,130.5)).mean('lon'), dpcentdt_sasia, 'kx', mew=2, ms=10, alpha=0.5)
ax2.plot(data.p_cent.sel(lon=np.arange(90.5,130.5)).mean('lon'), dpcentdt_sasia, alpha=0.3, color='k', linewidth=2)
#ax2.plot(data.p_cent.sel(lon=np.arange(110.5,120.5)).mean('lon'), dpcentdt_sasia, 'k')
ax2.plot(data_2p5.p_cent.mean('lon'), dpcentdt_2p5, linewidth=2)
ax2.plot(data_10.p_cent.mean('lon'), dpcentdt_10, linewidth=2)
ax2.plot(data_15.p_cent.mean('lon'), dpcentdt_15, linewidth=2)
ax2.plot(data_20.p_cent.mean('lon'), dpcentdt_20, linewidth=2)
ax2.set_title('South China Sea (90-130E)')

ax3.plot(data.p_cent.sel(lon=np.arange(-70.5,-45.5)).mean('lon'), dpcentdt_sa, 'kx', mew=2, ms=10, alpha=0.5)
ax3.plot(data.p_cent.sel(lon=np.arange(-70.5,-45.5)).mean('lon'), dpcentdt_sa, alpha=0.3, color='k', linewidth=2)
ax3.plot(data_2p5.p_cent.mean('lon'), dpcentdt_2p5, linewidth=2)
ax3.plot(data_10.p_cent.mean('lon'), dpcentdt_10, linewidth=2)
ax3.plot(data_15.p_cent.mean('lon'), dpcentdt_15, linewidth=2)
ax3.plot(data_20.p_cent.mean('lon'), dpcentdt_20, linewidth=2)
ax3.set_title('South America (45-70W)')

ax4.plot(data.p_cent.sel(lon=np.arange(17.5,35.5)).mean('lon'), dpcentdt_africa, 'kx', mew=2, ms=10, alpha=0.5)
ax4.plot(data.p_cent.sel(lon=np.arange(17.5,35.5)).mean('lon'), dpcentdt_africa, alpha=0.3, color='k', linewidth=2)
ax4.plot(data_2p5.p_cent.mean('lon'), dpcentdt_2p5, linewidth=2)
ax4.plot(data_10.p_cent.mean('lon'), dpcentdt_10, linewidth=2)
ax4.plot(data_15.p_cent.mean('lon'), dpcentdt_15, linewidth=2)
ax4.plot(data_20.p_cent.mean('lon'), dpcentdt_20, linewidth=2)
ax4.set_title('Africa (17-35E)')

ax3.set_xlabel('ITCZ latitude')
ax4.set_xlabel('ITCZ latitude')
ax1.set_ylabel('ITCZ migration rate')
ax3.set_ylabel('ITCZ migration rate')
    
for ax in [ax1,ax2,ax3,ax4]:
    ax.grid(True,linestyle=':')

plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)

# Save as a pdf
plt.savefig(plot_dir + 'bowties_gpcp.pdf', format='pdf')
plt.close()


# Set figure parameters
rcParams['figure.figsize'] = 7, 5
rcParams['font.size'] = 14
data.p_cent.sel(lon=np.arange(70.5,90.5)).mean('lon').plot(linewidth=2)
data.p_cent.sel(lon=np.arange(90.5,130.5)).mean('lon').plot(linewidth=2)
data.p_cent.sel(lon=np.arange(-70.5,-45.5)).mean('lon').plot(linewidth=2)
data.p_cent.sel(lon=np.arange(17.5,35.5)).mean('lon').plot(linewidth=2)
plt.xlabel('Pentad')
plt.ylabel('ITCZ latitude')
plt.legend(['India', 'South China Sea', 'South America','Africa'], loc='lower center')
plt.grid(True,linestyle=':')
plt.ylim(-20,20)
plt.xlim(0,73)
# Save as a pdf
plt.savefig(plot_dir + 'gpcp_itcz_lats.pdf', format='pdf')
plt.close()


rcParams['figure.figsize'] = 14, 6
rcParams['font.size'] = 14

fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4)
ax_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
j=0

for i in range(20,36,2):
    ax = ax_list[j]
    data.precipitation.sel(xofyear=i).plot.contourf(x='lon', y='lat', ax=ax, extend='max', cmap='Blues', levels = np.arange(2.,15.,2.), add_labels=False, add_colorbar=False)
    land_mask.plot.contour(x='lon', y='lat', ax=ax, levels=np.arange(-1.,2.,1.), add_labels=False, colors='0.7')
    ax.set_ylim([-30,30])
    ax.set_xlim([70,130])
    ax.set_yticks(np.arange(-30.,30.,10.))
    data.p_cent.sel(xofyear=i).plot.line(ax=ax, color='k', linewidth=2)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('')
    ax.grid(True,linestyle=':')
    j = j+1

ax1.set_ylabel('Latitude')
ax5.set_ylabel('Latitude')

for ax in ax_list[4:]:
    ax.set_xlabel('Longitude')

plt.savefig(plot_dir + 'gpcp_asia_onset.pdf', format='pdf')
plt.close()
    