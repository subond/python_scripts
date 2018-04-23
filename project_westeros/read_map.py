'''Load up a png of the Game of Thrones map and use to generate a land and topography mask'''

import scipy.misc as scp
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from netCDF4 import Dataset

# Set some constants
RADIUS = 6371.e3  # Assume radius is the same as Earth
#CIRCUMFERENCE = 2*np.pi*radius  #sanity check of distances vs circumference of planet
MILE2METRE = 1.60934*1000. # Conversion factor for miles to metres

# Import low resolution png of Game of Thrones map
image_in = scp.imread('got_continents_lowres_edits.png',flatten=True)

# Create mask of land and ocean using colour values
got_land = np.zeros((540,708))
got_land[image_in[::-1,:]>23.] = 1.

#Ratio of scale bar = 1200miles to total size of picture (width, height) is 2.1:(15,10)cm
#In this aspect wall is 1.5cm from top of picture
scalebar_m = 1200 * MILE2METRE  # Get scalebar length in metres
pic_width_m = scalebar_m / 2.1 * 15   # Width of whole picture in metres
pic_height_m = scalebar_m / 2.1 * 10  # Height of whole picture in metres

# Convert to lat-lon to get fraction of globe map covers
pic_width_lon = (pic_width_m / RADIUS) * (180./np.pi)
pic_height_lat = (pic_height_m / RADIUS) * (180./np.pi)

# Divide lat-lon range of map by shape to get the lat-lon size of each pixel
pixel_lon = pic_width_lon / got_land.shape[1]
pixel_lat = pic_height_lat / got_land.shape[0]

#Create lat lon axes corresponding to pixels
#Assume wall to be at 60N cf http://ibbenesecartographer.blogspot.co.uk/2013/05/the-true-size-of-north.html
wall_m_from_top = scalebar_m/2.1*1.5
wall_deg_from_top = (wall_m_from_top/RADIUS) * (180./np.pi)

lat_at_top = wall_deg_from_top + 60.
lat_at_bottom = lat_at_top - pic_height_lat

gotlat = np.arange(lat_at_bottom+pixel_lat, lat_at_top+pixel_lat, pixel_lat)
gotlon = np.arange(0., pic_width_lon, pixel_lon)

# We can now define a local landmask with lat-lon coords
got_land = xr.DataArray( got_land, [('gotlat', gotlat ), ('gotlon', gotlon )])

# Open up an existing land file to get model resolution
land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')


def find_equiv(x,y):
    '''For two coordinate axes x and y find the index in y representing the closest value to x'''
    x1=np.zeros(len(x))
    for i in range(0,len(x)):
        x1[i] = np.argmin(abs(y-x[i])) 
    x1=x1.astype(int)
    return x1

# Match GoT lat-lons to T42 lat-lons
lons_equiv = find_equiv(land.lon.values, gotlon)
lats_equiv = find_equiv(land.lat.values, gotlat)

# T42 GoT land mask in xarray
got_land_mask = got_land[lats_equiv, lons_equiv].values
got_land_mask = xr.DataArray( got_land_mask, land.coords)

# Mask off the bits that were outside the picture
lons = got_land_mask.lon[got_land_mask.lon >= max(gotlon)]
lats = got_land_mask.lat[ np.array(got_land_mask.lat >= max(gotlat)) | np.array(got_land_mask.lat < min(gotlat)) ]
got_land_mask.loc[dict(lon = lons)] = 0.
got_land_mask.loc[dict(lat = lats)] = 0.

# Plot to allow checking
got_land_mask.plot.contourf()
plt.plot([0,30],[60.,60.],'k')  # Indication of wall location
plt.plot([0,360.],[0.,0.],'w')  # Equator
plt.show()


#plt.contourf(got_land)
#plt.plot([0,707],[81,81],'k')  # Indication of wall location
#plt.show()


#Write land array to file
topo_filename = 'westeros.nc'
topo_file = Dataset(topo_filename, 'w', format='NETCDF3_CLASSIC')
nlon = got_land_mask.lon.shape[0]
nlat = got_land_mask.lat.shape[0]  
lat = topo_file.createDimension('lat', nlat)
lon = topo_file.createDimension('lon', nlon)
latitudes = topo_file.createVariable('lat','f4',('lat',))
longitudes = topo_file.createVariable('lon','f4',('lon',))
land_array_netcdf = topo_file.createVariable('land_mask','f4',('lat','lon',))
latitudes[:] = got_land_mask.lat.values
longitudes[:] = got_land_mask.lon.values
land_array_netcdf[:] = got_land_mask.values
topo_file.close()
print 'Output written to: ' + topo_filename

