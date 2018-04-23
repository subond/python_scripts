# Evaluate the mass streamfunction as a function of time, look at progression of peak throughout year

from physics import mass_streamfunction
from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

data = time_means('sn_3.000', [145,469], filename='atmos_pentad', timeav='pentad',period_fac=3.)
psi_30_3 = mass_streamfunction(data, a=6376.0e3)/1e9
print 'psi evaluated'

data_ap10 = time_means('aquaplanet_10m', [109,145], filename='atmos_pentad', timeav='pentad',period_fac=1.)
psi_10_1 = mass_streamfunction(data_ap10, a=6376.0e3)/1e9

#u = data.ucomp
#print 'u loaded'

psi_max_30_3 = np.argmin( psi_30_3.values[0:40,:,:], axis=0)
psi_max_10_1 = np.argmin( psi_10_1.values[0:40,:,:], axis=0)
#36
#psi_30_3[:,:,36].plot.contourf(x='xofyear', y='lat', levels = np.arange(-450,451,150), extend='neither')

plt.plot(data.xofyear/3., data.lat[psi_max_30_3[:,27]],'k')
plt.plot(data_ap10.xofyear, data.lat[psi_max_10_1[:,27]],'r')
plt.ylim(-45,45)
plt.xlabel('Pentad (normalised to Earth)')
plt.ylabel('Latitude')
plt.legend(['3x','1x'])
#plt.xlim(64,136)

#plt.savefig('/scratch/rg419/plots/seasons_and_rotation/psi_3times_30m.png')
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/psi_110_330.png')
plt.close()


plt.plot(data_ap10.xofyear, data.lat[psi_max_30_3[61:133,27]],'k')
plt.plot(data_ap10.xofyear, data.lat[psi_max_10_1[:,27]],'r')
plt.ylim(-45,45)
plt.xlabel('Pentad')
plt.ylabel('Latitude')
plt.legend(['3x','1x'])
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/psi_110_330_zoom.png')
plt.close()