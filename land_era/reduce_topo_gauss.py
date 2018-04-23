# 08/02/2018 Experimenting with using the Gaussian topography/it's mask to remove Tibet from experiments. Conclusion - setting topography to 0 seems odd, the rest of the continent is elevated. For now just go with plan a and set it to 500.

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
resolution_file = Dataset(GFDL_BASE + 'src/extra/python/scripts/gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')
lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]
lonb = resolution_file.variables['lonb'][:]
latb = resolution_file.variables['latb'][:]
nlon=lons.shape[0]
nlat=lats.shape[0]
topo_array = np.zeros((nlat,nlon))
land_array = np.zeros((nlat,nlon))

#make 2d arrays of latitude and longitude
lon_array, lat_array = np.meshgrid(lons,lats)
lonb_array, latb_array = np.meshgrid(lonb,latb)
    
    
#Tibet from Sauliere 2012
h_0 = 1.
central_lon = 82.5
central_lat = 28
L_1 = 12.5
L_2 = 12.5
gamma_1 = -49.5
gamma_2 = -15. #-18.
delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/L_1
delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/L_2
h_arr_tibet_no_amp = np.exp(-(delta_1**2.))*(1./delta_2)*np.exp(-0.5*(np.log(delta_2))**2.)
maxval = np.nanmax(h_arr_tibet_no_amp) #For some reason my maximum value of h_arr_tibet_no_amp > 1. Renormalise so h_0 sets amplitude. 
h_arr_tibet = (h_arr_tibet_no_amp/maxval)*h_0
idx_tibet = (h_arr_tibet / h_0 > 0.05)

topo_array = np.zeros((nlat,nlon))
topo_array[idx_tibet] =  h_arr_tibet[idx_tibet]
topo_array = 1 - idx_tibet

data = xr.open_dataset('/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
data.load()

topo_out = data.zsurf * topo_array
land_out = data.land_mask

#Write land and topography arrays to file
topo_file = Dataset('reduce_topo_gauss_test5.nc', 'w', format='NETCDF3_CLASSIC')
lat = topo_file.createDimension('lat', nlat)
lon = topo_file.createDimension('lon', nlon)
latitudes = topo_file.createVariable('lat','f4',('lat',))
longitudes = topo_file.createVariable('lon','f4',('lon',))
topo_array_netcdf = topo_file.createVariable('zsurf','f4',('lat','lon',))
land_array_netcdf = topo_file.createVariable('land_mask','f4',('lat','lon',))
latitudes[:] = lats
longitudes[:] = lons
topo_array_netcdf[:] = topo_out.values
land_array_netcdf[:] = land_out.values
topo_file.close()