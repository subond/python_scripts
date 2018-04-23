# Evaluate angular momentum 

from physics import mass_streamfunction
from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

a=6376.0e3
omega = 7.2921150e-5

data = time_means('aquaplanet_2m', [121,145], filename='atmos_daily', timeav='pentad')
u = data.ucomp.mean('lon')
print 'u loaded'

a_cos = a*np.cos(data.lat * np.pi/180.)
m = (u + omega*a_cos)*a_cos

m[35,:,:].plot.contourf(x='lat',y='pfull', yincrease=False)
plt.xlim(-45,45)
plt.figure(2)
m[25,:,:].plot.contourf(x='lat',y='pfull', yincrease=False)
plt.xlim(-45,45)
plt.show()