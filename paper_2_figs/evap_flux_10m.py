"""Plot the evaporative flux"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from climatology import precip_centroid

    
rcParams['figure.figsize'] = 10, 5
rcParams['font.size'] = 20
    
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')

# Get total precip
try:
    data['precipitation'] = data.condensation_rain + data.convection_rain
except:
    data['precipitation'] = data.precipitation
    
precip_temp = np.zeros(data.precipitation.values.shape)
lhe_temp = np.zeros(data.flux_lhe.values.shape)

n = len(data.xofyear.values)/2

for i in range(0,n):
    precip_temp[i,:,:] = (data.precipitation[i,:,:].values + data.precipitation[i+n,::-1,:].values)/2.
    precip_temp[i+n,:,:] = precip_temp[i,::-1,:]
    lhe_temp[i,:,:] = (data.flux_lhe[i,:,:].values + data.flux_lhe[i+n,::-1,:].values)/2.
    lhe_temp[i+n,:,:] = lhe_temp[i,::-1,:]
    
precip_temp = xr.DataArray(precip_temp, coords=[data.xofyear.values, data.lat, data.lon], dims=['xofyear', 'lat', 'lon'])
lhe_temp = xr.DataArray(lhe_temp, coords=[data.xofyear.values, data.lat, data.lon], dims=['xofyear', 'lat', 'lon'])

data['precipitation'] = precip_temp
data['flux_lhe'] = lhe_temp


# Locate precipitation centroid
precip_centroid(data)

    
data = data.mean('lon')
levels_flux = np.arange(-300.,305.,20.)

fig = plt.figure()
ax1 = plt.subplot(111)

f1=(-1.*data.flux_lhe).plot.contourf(x='xofyear', y='lat', ax=ax1, levels=levels_flux, extend = 'both', add_labels=False, cmap='RdBu_r', add_colorbar=False)
data.p_cent.plot(ax=ax1, color='k', linewidth=2) 

ax1.set_ylim([-60.,60.])
ax1.set_ylabel('Latitude')
ax1.set_xlabel('Pentad')
ax1.grid(True,linestyle=':')

cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'vertical',fraction=0.15, pad=0.05, aspect=30)
cb1.set_label('Evap. heat flux, W/m^2')

plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.18)

plt.savefig(plot_dir + 'evap_flux_10m.pdf', format='pdf')
plt.close()

