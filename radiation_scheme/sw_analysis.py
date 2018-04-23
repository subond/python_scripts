# Compare orientation of AM and Psi surfaces over year

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np


name_temp = '/scratch/rg419/Data_moist/rad_scheme/run%03d/atmos_daily.nc'
names = [name_temp % m for m in range( 109, 481)  ]
#read data into xarray 
#data = xr.open_mfdataset( names, decode_times=False, chunks={'time': 30})
data = xr.open_dataset( '/scratch/rg419/Data_moist/rad_scheme/run013/atmos_daily.nc', decode_times=False)


dtau = -np.log(data.swdflx.values[:,1:41, ...]/data.swdflx.values[:,0:40, ...])

data['dtau'] =  (('time','pfull','lat','lon'), dtau)	

plt.plot(data.sphum[0,10,...].values.flatten(), data.dtau[0,10,...].values.flatten(), 'x')
plt.ylim([0,0.1])
plt.show()