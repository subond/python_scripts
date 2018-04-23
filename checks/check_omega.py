#import packages
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


#open files
nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/cssp_run_topo_reg/np16/run1/atmos_daily.nc'
fh = Dataset(nc_file, mode='r')
lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
pfulls = fh.variables['pfull'][:]
times = fh.variables['time'][:]
w_j = fh.variables['omega'][:]
fh.close()
nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/cssp_run_topo_reg/np16/run2/atmos_daily.nc'
fh = Dataset(nc_file, mode='r')
w_f = fh.variables['omega'][:]
fh.close()

land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/cssp_run_topo_reg/input/land.nc'
lfh = Dataset(land_file, 'r', format='NETCDF3_CLASSIC')
land_mask = lfh.variables['land_mask'][:]
lfh.close()

# Get some parameters for the Stereographic Projection
lon_0 = lons.mean()
lat_0 = lats.mean()

lon, lat = np.meshgrid(lons, lats)
m = Basemap(lat_0=lat_0,lon_0=lon_0)
xi, yi = m(lon, lat)

plt.figure(1)
plt.set_cmap('bwr')
cs = m.contourf(xi,yi,np.squeeze(np.mean( (w_j[:,37,:,:] + w_f[:,37,:,:])/2. ,0)), np.arange(-0.2,0.21,0.01))
plt.contour(xi,yi,land_mask,[0,1],colors='k')
# Add Grid Lines
m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")


#open files
nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/cssp_run_topo/np16/run1/atmos_daily.nc'
fh = Dataset(nc_file, mode='r')
lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
pfulls = fh.variables['pfull'][:]
times = fh.variables['time'][:]
w_j = fh.variables['omega'][:]
fh.close()
nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/cssp_run_topo/np16/run2/atmos_daily.nc'
fh = Dataset(nc_file, mode='r')
w_f = fh.variables['omega'][:]
fh.close()

land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/cssp_run_topo/input/land.nc'
lfh = Dataset(land_file, 'r', format='NETCDF3_CLASSIC')
land_mask = lfh.variables['land_mask'][:]
lfh.close()

# Get some parameters for the Stereographic Projection
lon_0 = lons.mean()
lat_0 = lats.mean()

lon, lat = np.meshgrid(lons, lats)
m = Basemap(lat_0=lat_0,lon_0=lon_0)
xi, yi = m(lon, lat)

plt.figure(2)
cs = m.contourf(xi,yi,np.squeeze(np.mean( (w_j[:,37,:,:] + w_f[:,37,:,:])/2. ,0)), np.arange(-0.2,0.21,0.01))
plt.contour(xi,yi,land_mask,[0,1],colors='k')
# Add Grid Lines
m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
# Add Colorbar
cbar = m.colorbar(cs, location='bottom')


plt.show()


