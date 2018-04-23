# 26/01/2018 Load in ERA T42 land mask and reduce topo height over selected regions.

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import scipy.interpolate as scip
import xarray as xr


def reduce_topo(filename, reduce_parts=['None'], check=False, thresh=0.):

    data = xr.open_dataset('/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    data.load()

    land_out = data.land_mask.copy()
    topo_out = data.zsurf.copy()
    
    if 'america' in reduce_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= 190 and data.lon[i] < 347]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= -60]
        topo_section = topo_out[min(lats):max(lats), min(lons):max(lons)].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats):max(lats), min(lons):max(lons)] = topo_section
    
    if 'south_america' in reduce_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= 280 and data.lon[i] < 330]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= -60 and data.lat[i] < 15]
        topo_section = topo_out[min(lats):max(lats), min(lons):max(lons)].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats):max(lats), min(lons):max(lons)] = topo_section
    
    if 'north_america' in reduce_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= 190 and data.lon[i] < 347]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= 15]
        lons_b = [i for i in range(len(data.lon)) if data.lon[i] >= 270 and data.lon[i] < 280]
        lats_b = [i for i in range(len(data.lat)) if data.lat[i] >= 6]
        topo_section = topo_out[min(lats):max(lats), min(lons):max(lons)].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats):max(lats), min(lons):max(lons)] = topo_section
        
        topo_section = topo_out[min(lats_b):max(lats_b), min(lons_b):max(lons_b)].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats_b):max(lats_b), min(lons_b):max(lons_b)] = topo_section

        
    if 'australia' in reduce_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= 95 and data.lon[i] < 200]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= -60 and data.lat[i] < 7]
        topo_section = topo_out[min(lats):max(lats), min(lons):max(lons)].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats):max(lats), min(lons):max(lons)] = topo_section
    
    if 'africa' in reduce_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] < 65]
        lons_b = [i for i in range(len(data.lon)) if data.lon[i] >= 345]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= -50 and data.lat[i] < 27]
        topo_section = topo_out[min(lats):max(lats), min(lons):max(lons)].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats):max(lats), min(lons):max(lons)] = topo_section
        
        topo_section = topo_out[min(lats):max(lats), min(lons_b):].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats):max(lats), min(lons_b):] = topo_section
        
    
    if 's_asia' in reduce_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= 65 and data.lon[i] < 130]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= 6 and data.lat[i] < 27]
        topo_section = topo_out[min(lats):max(lats), min(lons):max(lons)].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats):max(lats), min(lons):max(lons)] = topo_section
    
    if 'europe' in reduce_parts:
        lons = [i for i in range(len(data.lon)) if data.lon[i] < 193]
        lons_b = [i for i in range(len(data.lon)) if data.lon[i] >= 345]
        lats = [i for i in range(len(data.lat)) if data.lat[i] >= 26]
        topo_section = topo_out[min(lats):max(lats), min(lons):max(lons)].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats):max(lats), min(lons):max(lons)] = topo_section
        
        topo_section = topo_out[min(lats):max(lats), min(lons_b):].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats):max(lats), min(lons_b):] = topo_section
        
    if 'antarctica' in reduce_parts:
        lats = [i for i in range(len(data.lat)) if data.lat[i] < -60]
        topo_section = topo_out[min(lats):max(lats), min(lons):max(lons)].values
        topo_section[topo_section>thresh] = thresh
        topo_out[min(lats):max(lats), min(lons):max(lons)] = topo_section
        
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
    lonb = resolution_file.variables['lonb'][:]
    latb = resolution_file.variables['latb'][:]
    nlon=lons.shape[0]
    nlat=lats.shape[0]    
    topo_array = np.zeros((nlat,nlon))
    land_array = np.zeros((nlat,nlon))
    
    #Write land and topography arrays to file
    topo_file = Dataset(filename, 'w', format='NETCDF3_CLASSIC')
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
    print 'Output written to: ' + filename
    



reduce_topo('land_era_no_TIP_0.nc', reduce_parts=['europe', 's_asia'])




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



