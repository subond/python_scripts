#load in topography output from model initialisation and compare with input from land.nc

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#set run name
run_fol = 'monsoon_unsmooth/np16'
inp_fol = 'monsoon_unsmooth'

#set run name
#run_fol = 'cssp_run_topo/np16'
#inp_fol = 'cssp_run_topo'

#open files

nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/' + run_fol + '/run1/atmos_daily.nc'

fh = Dataset(nc_file, mode='r')
lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
p = fh.variables['pfull'][:]
times = fh.variables['time'][:]
zsurf = fh.variables['zsurf'][:]
fh.close()

land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/' + inp_fol + '/input/land.nc'
lfh = Dataset(land_file, 'r', format='NETCDF3_CLASSIC')
land_array = lfh.variables['land_mask'][:]
zsurf_in = lfh.variables['zsurf'][:]
lfh.close()

lon_array, lat_array = np.meshgrid(lons,lats)
lon_0 = lons.mean()
lat_0 = lats.mean()
m = Basemap(lat_0=lat_0,lon_0=lon_0)
xi, yi = m(lon_array, lat_array)

plt.figure(1)
m.contour(xi,yi,land_array)
cs = m.contourf(xi,yi,zsurf_in,np.arange(-600,6500,600), cmap=plt.get_cmap('RdBu_r'))
m.drawparallels(np.arange(-60., 61., 30.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)
m.drawcoastlines()
cb = plt.colorbar(cs, shrink=0.5, extend='both')
plt.savefig('zsurf_in_us.png')

plt.figure(2)
m.contour(xi,yi,land_array)
cs = m.contourf(xi,yi,zsurf,np.arange(-600,6500,600), cmap=plt.get_cmap('RdBu_r'))
m.drawparallels(np.arange(-60., 61., 30.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)
m.drawcoastlines()
cb = plt.colorbar(cs, shrink=0.5, extend='both')
plt.savefig('zsurf_us.png')

plt.figure(3)
cs = m.contourf(xi,yi,zsurf-zsurf_in,np.arange(-2000,2000,200), cmap=plt.get_cmap('RdBu_r'))
m.drawparallels(np.arange(-60., 61., 30.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)
m.drawcoastlines()
cb = plt.colorbar(cs, shrink=0.5, extend='both')
plt.savefig('zsurf_outmin_us.png')

plt.show()