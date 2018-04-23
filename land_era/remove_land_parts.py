# 18/01/2018 Load in ERA T42 land mask and remove selected parts.

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import scipy.interpolate as scip
import xarray as xr


def remove_land_parts(filename, var_name='land_mask', add_topo=True, remove_parts=['None'], check=False, use_boundaries=False):

    data = xr.open_dataset('/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    data.load()

    land_out = data.land_mask.copy()
    topo_out = data.zsurf.copy()
    
    if 'america' in remove_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= 190 and data.lon[i] < 347]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= -60]
        land_out[min(lats):max(lats), min(lons):max(lons)] = 0.
        topo_out[min(lats):max(lats), min(lons):max(lons)] = 0.
    
    if 'south_america' in remove_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= 280 and data.lon[i] < 330]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= -60 and data.lat[i] < 15]
        land_out[min(lats):max(lats), min(lons):max(lons)] = 0.
        topo_out[min(lats):max(lats), min(lons):max(lons)] = 0.
    
    if 'north_america' in remove_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= 190 and data.lon[i] < 347]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= 15]
        lons_b = [i for i in range(len(data.lon)) if data.lon[i] >= 270 and data.lon[i] < 280]
        lats_b = [i for i in range(len(data.lat)) if data.lat[i] >= 6]
        land_out[min(lats):max(lats), min(lons):max(lons)] = 0.
        land_out[min(lats_b):max(lats_b), min(lons_b):max(lons_b)] = 0.
        topo_out[min(lats):max(lats), min(lons):max(lons)] = 0.
        topo_out[min(lats_b):max(lats_b), min(lons_b):max(lons_b)] = 0.
        
    if 'australia' in remove_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= 95 and data.lon[i] < 200]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= -60 and data.lat[i] < 7]
        land_out[min(lats):max(lats), min(lons):max(lons)] = 0.
        topo_out[min(lats):max(lats), min(lons):max(lons)] = 0.
    
    if 'africa' in remove_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] < 65]
        lons_b = [i for i in range(len(data.lon)) if data.lon[i] >= 345]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= -50 and data.lat[i] < 27]
        land_out[min(lats):max(lats), min(lons):max(lons)] = 0.
        land_out[min(lats):max(lats), min(lons_b):] = 0.
        topo_out[min(lats):max(lats), min(lons):max(lons)] = 0.
        topo_out[min(lats):max(lats), min(lons_b):] = 0.
    
    if 's_asia' in remove_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= 65 and data.lon[i] < 130]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= 6 and data.lat[i] < 27]
        land_out[min(lats):max(lats), min(lons):max(lons)] = 0.
        topo_out[min(lats):max(lats), min(lons):max(lons)] = 0.
    
    if 'europe' in remove_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] < 193]
        lons_b = [i for i in range(len(data.lon)) if data.lon[i] >= 345]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= 26]
        land_out[min(lats):max(lats), min(lons):max(lons)] = 0.
        land_out[min(lats):max(lats), min(lons_b):] = 0.
        topo_out[min(lats):max(lats), min(lons):max(lons)] = 0.
        topo_out[min(lats):max(lats), min(lons_b):] = 0.
        
    if 'antarctica' in remove_parts:
        lats = [i for i in range(len(data.lat)) if data.lat[i] < -60]
        land_out[min(lats):max(lats), :] = 0.
        topo_out[min(lats):max(lats), :] = 0.
        
    if check:
        plt.figure(1)
        data.zsurf.plot.contourf(levels=np.arange(0,6000,250))
        data.land_mask.plot.contour()
    
        plt.figure(2)
        topo_out.plot.contourf(levels=np.arange(0,6000,250))
        land_out.plot.contour()
    
        plt.figure(3)
        (data.zsurf - topo_out).plot.contourf(levels=np.arange(0,6000,250))
        (data.land_mask - land_out).plot.contour()
        
        plt.show()
    
    # Common features of set-ups
    # specify resolution
    t_res = 42
    #read in grid from approriate file
    GFDL_BASE = os.environ['GFDL_BASE']
    resolution_file = Dataset(GFDL_BASE + 'src/extra/python/scripts/gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')
    lons = resolution_file.variables['lon'][:]
    lats = resolution_file.variables['lat'][:]
    lonbs = resolution_file.variables['lonb'][:]
    latbs = resolution_file.variables['latb'][:]
    nlon=lons.shape[0]
    nlat=lats.shape[0]   
    nlonb=lonbs.shape[0]
    nlatb=latbs.shape[0] 
    land_array = np.zeros((nlat,nlon))
    
    #Write land and topography arrays to file
    topo_file = Dataset(filename, 'w', format='NETCDF3_CLASSIC')
    lat = topo_file.createDimension('lat', nlat)
    lon = topo_file.createDimension('lon', nlon)
    latitudes = topo_file.createVariable('lat','f4',('lat',))
    longitudes = topo_file.createVariable('lon','f4',('lon',))


    if use_boundaries:
        latb = topo_file.createDimension('latb', nlatb)
        latitudebs = topo_file.createVariable('latb','d',('latb',))
        latitudes.edges = 'latb'
        latitudebs.units = 'degrees_N'.encode('utf-8')
        latitudebs.cartesian_axis = 'Y'
        latitudebs.long_name = 'latitude edges'
        latitudebs[:] = latbs
        latitudes.units = 'degrees_N'.encode('utf-8')
        latitudes.cartesian_axis = 'Y'
        latitudes.long_name = 'latitude'
        
        lonb = topo_file.createDimension('lonb', nlonb)
        longitudebs = topo_file.createVariable('lonb','d',('lonb',))
        longitudes.edges = 'lonb'
        longitudebs.units = 'degrees_E'.encode('utf-8')
        longitudebs.cartesian_axis = 'X'
        longitudebs.long_name = 'longitude edges'
        longitudebs[:] = lonbs
        longitudes.units = 'degrees_E'.encode('utf-8')
        longitudes.cartesian_axis = 'X'
        longitudes.long_name = 'longitude'
    
    if add_topo:
        topo_array = np.zeros((nlat,nlon))    
        topo_array_netcdf = topo_file.createVariable('zsurf','f4',('lat','lon',))
        topo_array_netcdf[:] = topo_out.values
    
    land_array_netcdf = topo_file.createVariable(var_name,'f4',('lat','lon',))
    latitudes[:] = lats
    longitudes[:] = lons
    land_array_netcdf[:] = land_out.values
    topo_file.close()
    print 'Output written to: ' + filename
    



remove_land_parts('icy_america.nc', var_name='icy_america', add_topo=False, remove_parts=['europe', 'antarctica', 'africa', 's_asia', 'australia'], use_boundaries=True)

remove_land_parts('test.nc')




#land_array = era.lsm[0,lats_equiv,lons_equiv].values
#topo_array = era.z[0,lats_equiv,lons_equiv].values/9.8
#print land_array
#idx = (land_array == 0.) & (topo_array != 0.)
#topo_array[idx] = 0. 
        
#print land_array.shape,topo_array.shape
#plt.contourf(topo_array)
#plt.contour(land_array)
#plt.show()

#Write land and topography arrays to file
#topo_filename = 'land.nc'
#topo_file = Dataset(topo_filename, 'w', format='NETCDF3_CLASSIC')
#lat = topo_file.createDimension('lat', nlat)
#lon = topo_file.createDimension('lon', nlon)
#latitudes = topo_file.createVariable('lat','f4',('lat',))
#longitudes = topo_file.createVariable('lon','f4',('lon',))
#topo_array_netcdf = topo_file.createVariable('zsurf','f4',('lat','lon',))
#land_array_netcdf = topo_file.createVariable('land_mask','f4',('lat','lon',))
#latitudes[:] = lats
#longitudes[:] = lons
#topo_array_netcdf[:] = topo_array
#land_array_netcdf[:] = land_array
#topo_file.close()
#print 'Output written to: ' + topo_filename



