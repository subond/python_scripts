# Compare orientation of AM and Psi surfaces over year

from physics import mass_streamfunction
from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from finite_difference import cfd


a=6376.0e3
omega = 7.2921150e-5

data = time_means('topo_10m', [121,145], filename='atmos_daily', timeav='pentad')
u = data.ucomp.mean('lon')
print 'u loaded'

a_cos = a*np.cos(data.lat * np.pi/180.)
m = (u + omega*a_cos)*a_cos
psi = abs(mass_streamfunction(data, a=6376.0e3)/1e9)


dpsidy = xr.DataArray( cfd( psi.values, psi.lat*np.pi/180,0),   [('lat', psi.lat), ('xofyear', psi.xofyear)  , ('pfull', psi.pfull )])
dpsidp = xr.DataArray( cfd( psi.values, psi.pfull*100,2),    [('lat', psi.lat), ('xofyear', psi.xofyear)  , ('pfull', psi.pfull )])

dmdy = xr.DataArray( cfd( m.values, m.lat*np.pi/180,2),   [('xofyear', psi.xofyear), ('pfull', m.pfull ), ('lat', m.lat)])
dmdp = xr.DataArray( cfd( m.values, m.pfull*100,1),   [('xofyear', psi.xofyear), ('pfull', m.pfull ), ('lat', m.lat)])

dpsi_mag = np.sqrt(dpsidy*dpsidy + dpsidp*dpsidp)
dm_mag = np.sqrt(dmdy*dmdy + dmdp*dmdp)

dpsidy = dpsidy/dpsi_mag
dpsidp = dpsidp/dpsi_mag
dmdy = dmdy/dm_mag
dmdp = dmdp/dm_mag

j_diff = dmdy - dpsidy
k_diff = dmdp - dpsidp

m_psi = (dmdy*dpsidy + dmdp*dpsidp)#/dm_mag/dpsi_mag

plt.figure(0)
psi[:,34,:].plot.contourf(x='lat', y='pfull', yincrease=False)
plt.xlim(-45,45)

plt.figure(1)
m[34,:,:].plot.contourf(x='lat', y='pfull', yincrease=False)
plt.xlim(-45,45)

plt.figure(2)
j_diff[:,27,:].plot.contourf(x='xofyear', y='lat', levels=np.arange(-0.5,0.55,0.05))
plt.ylim(-45,45)


plt.show()

