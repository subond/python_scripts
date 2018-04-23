import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import scipy.interpolate as scip
import xarray as xr


# Common features of set-ups
# specify resolution
t_res = 42
#read in grid from approriate file
GFDL_BASE = os.environ['GFDL_BASE']
resolution_file = Dataset('/scratch/rg419/Data_moist/ap_2_t85/run001/atmos_monthly.nc', 'r', format='NETCDF3_CLASSIC')
lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]
lonb = resolution_file.variables['lonb'][:]
latb = resolution_file.variables['latb'][:]
nlon=lons.shape[0]
nlat=lats.shape[0]    
topo_array = np.zeros((nlat,nlon))
land_array = np.zeros((nlat,nlon))

era = xr.open_dataset('ERA-I_Invariant_0125.nc')

def find_equiv(x,y):
    #for two coordinate axes x and y find the index in y representing the closest value to x
    x1=np.zeros(len(x))
    for i in range(0,len(x)):
        x1[i] = np.argmin(abs(y-x[i])) 
    x1=x1.astype(int)
    return x1

lons_equiv = find_equiv(lons,era.longitude.values)
lats_equiv = find_equiv(lats,era.latitude.values)

land_array = era.lsm[0,lats_equiv,lons_equiv].values
topo_array = era.z[0,lats_equiv,lons_equiv].values/9.8
print land_array
idx = (land_array == 0.) & (topo_array != 0.)
topo_array[idx] = 0. 
        
print land_array.shape,topo_array.shape
plt.contourf(topo_array)
plt.contour(land_array)
plt.show()

#Write land and topography arrays to file
topo_filename = 'era_land_t85.nc'
topo_file = Dataset(topo_filename, 'w', format='NETCDF3_CLASSIC')
lat = topo_file.createDimension('lat', nlat)
lon = topo_file.createDimension('lon', nlon)
latitudes = topo_file.createVariable('lat','f4',('lat',))
longitudes = topo_file.createVariable('lon','f4',('lon',))
topo_array_netcdf = topo_file.createVariable('zsurf','f4',('lat','lon',))
land_array_netcdf = topo_file.createVariable('land_mask','f4',('lat','lon',))
latitudes[:] = lats
longitudes[:] = lons
topo_array_netcdf[:] = topo_array
land_array_netcdf[:] = land_array
topo_file.close()
print 'Output written to: ' + topo_filename



