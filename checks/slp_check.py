#load in topography output from model initialisation and compare with input from land.nc

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#set run name
run_fol = 'monsoon_12m/np8'
inp_fol = 'monsoon_12m'

#set run name
#run_fol = 'cssp_run_final/np16'
#inp_fol = 'cssp_run_final'

#open files

nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/' + run_fol + '/run1/atmos_daily.nc'

fh = Dataset(nc_file, mode='r')
lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
p = fh.variables['pfull'][:]
times = fh.variables['time'][:]
ps = fh.variables['ps'][:]
slp = fh.variables['slp'][:]
fh.close()

land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/' + inp_fol + '/input/land.nc'
#land_file = '/scratch/rg419/Data_moist/runfile_archive/' + inp_fol + '/input/land.nc'
lfh = Dataset(land_file, 'r', format='NETCDF3_CLASSIC')
land_array = lfh.variables['land_mask'][:]
lfh.close()

lon_array, lat_array = np.meshgrid(lons,lats)
lon_0 = lons.mean()
lat_0 = lats.mean()
m = Basemap(lat_0=lat_0,lon_0=lon_0)
xi, yi = m(lon_array, lat_array)

plt.figure(1)
cs = m.contourf(xi,yi,ps[9,:,:], np.arange(56000,104000,2000), cmap=plt.get_cmap('RdBu_r'))
m.drawparallels(np.arange(-60., 61., 30.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)
m.drawcoastlines()
cb = plt.colorbar(cs, shrink=0.5, extend='both')
plt.savefig('ps9.png')

plt.figure(2)
cs = m.contourf(xi,yi,ps[19,:,:]-ps[9,:,:], cmap=plt.get_cmap('RdBu_r'))
m.drawparallels(np.arange(-60., 61., 30.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)
m.drawcoastlines()
cb = plt.colorbar(cs, shrink=0.5, extend='both')
plt.savefig('ps19.png')

plt.figure(3)
cs = m.contourf(xi,yi,slp[9,:,:], np.arange(56000,104000,2000), cmap=plt.get_cmap('RdBu_r'))
m.drawparallels(np.arange(-60., 61., 30.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)
m.drawcoastlines()
cb = plt.colorbar(cs, shrink=0.5, extend='both')
plt.savefig('slp9.png')

plt.figure(4)
cs = m.contourf(xi,yi,slp[19,:,:]-slp[9,:,:], cmap=plt.get_cmap('RdBu_r'))
m.drawparallels(np.arange(-60., 61., 30.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)
m.drawcoastlines()
cb = plt.colorbar(cs, shrink=0.5, extend='both')
plt.savefig('slp19.png')


plt.figure(5)
cs = m.contourf(xi,yi,ps[9,:,:]-slp[9,:,:], cmap=plt.get_cmap('RdBu_r'))
m.drawparallels(np.arange(-60., 61., 30.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)
m.drawcoastlines()
cb = plt.colorbar(cs, shrink=0.5, extend='both')
plt.savefig('ps_m_slp.png')


plt.show()