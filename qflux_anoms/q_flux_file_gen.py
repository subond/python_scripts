# Generate a steady state q flux input file, with anomaly against a swamp ocean.

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from data_handling import cell_area


month_dict={'jan':0, 'feb':1, 'mar':2, 'apr':3, 'may':4, 'jun':5, 'jul':6, 'aug':7, 'sep':8, 'oct':9, 'nov':10, 'dec':11}

#read in grid from approriate file
resolution_file = Dataset('/scratch/rg419/Data_moist/climatologies/ss_92.000.nc', 'r', format='NETCDF3_CLASSIC')

lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]

lonbs = resolution_file.variables['lonb'][:]
latbs = resolution_file.variables['latb'][:]

nlon=lons.shape[0]
nlat=lats.shape[0]

nlonb=lonbs.shape[0]
nlatb=latbs.shape[0]

area = cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')   # Area of grid cells



warmpool_loc_list=[
[0., 95.]#, [5., 95.], [10., 95.], [15., 95.], [20., 95.], [25., 95.]
]


for file_values in warmpool_loc_list:

    warmpool_array = np.zeros([nlat,nlon])

    warmpool_lat_centre=file_values[0]
    warmpool_lon_centre=file_values[1]

    print warmpool_lat_centre,warmpool_lon_centre

    warmpool_width = 7.5 #15.
    warmpool_width_lon = 7.5 #45.
    warmpool_amp = 200.


    for j in np.arange(len(lats)):
         lat = 0.5*(latbs[j+1] + latbs[j])
         lat = (lat-warmpool_lat_centre)/warmpool_width
    #     if( np.absolute(lat) <= 1.0 ):
         for i in np.arange(len(lons)):
              lon = 0.5*(lonbs[i+1] + lonbs[i])
              lon = (lon-warmpool_lon_centre)/warmpool_width_lon
              if( lat**2.+lon**2. <= 1.0 ):
                  warmpool_array[j,i] = (1.-lat**2.-lon**2.)*warmpool_amp


    warmpool_integral = np.sum(area*warmpool_array)

    lon_array, lat_array = np.meshgrid(lons,lats)


    cooling_idx = (warmpool_array == 0.)

    cooling_array = np.zeros((nlat,nlon))
    cooling_array[cooling_idx] = 1.

    warmpool_array[cooling_idx]=-1.*warmpool_integral/np.sum(area[cooling_idx])

    warmpool_integral_adj = np.sum(area*warmpool_array)

    print 'warmpool_integral', warmpool_integral
    print 'warmpool_integral_adj', warmpool_integral_adj


    plt.figure()

    lon_0 = lons.mean()
    lat_0 = lats.mean()

    lon_array, lat_array = np.meshgrid(lons,lats)
    lonb_array, latb_array = np.meshgrid(lonbs,latbs)

    m = Basemap(lat_0=lat_0,lon_0=lon_0)
    xi, yi = m(lon_array, lat_array)

    cs = m.contourf(xi,yi, warmpool_array, cmap=plt.get_cmap('RdBu_r'))
    plt.xticks(np.linspace(0,360,13))
    plt.yticks(np.linspace(-90,90,7))
    cb = plt.colorbar(cs, shrink=0.5, extend='both')
    
    #plt.show()

    description_string = str(int(warmpool_lon_centre))+'_lat_'+str(int(warmpool_lat_centre))+'_a_'+str(int(warmpool_amp)) +'_small'
    
    var_name_out = 'qflux_' + description_string

    warmpool_file = Dataset(var_name_out+'.nc', 'w', format='NETCDF3_CLASSIC')
    lat = warmpool_file.createDimension('lat', nlat)
    lon = warmpool_file.createDimension('lon', nlon)

    latb = warmpool_file.createDimension('latb', nlatb)
    lonb = warmpool_file.createDimension('lonb', nlonb)

    latitudes = warmpool_file.createVariable('lat','f4',('lat',))
    longitudes = warmpool_file.createVariable('lon','f4',('lon',))

    latitudebs = warmpool_file.createVariable('latb','f4',('latb',))
    longitudebs = warmpool_file.createVariable('lonb','f4',('lonb',))

    latitudes.units = 'degrees_N'.encode('utf-8')
    latitudes.cartesian_axis = 'Y'
    latitudes.edges = 'latb'
    latitudes.long_name = 'latitude'

    longitudes.units = 'degrees_E'.encode('utf-8')
    longitudes.cartesian_axis = 'X'
    longitudes.edges = 'lonb'
    longitudes.long_name = 'longitude'

    latitudebs.units = 'degrees_N'.encode('utf-8')
    latitudebs.cartesian_axis = 'Y'
    latitudebs.long_name = 'latitude edges'

    longitudebs.units = 'degrees_E'.encode('utf-8')
    longitudebs.cartesian_axis = 'X'
    longitudebs.long_name = 'longitude edges'

    warmpool_array_netcdf = warmpool_file.createVariable('ocean_qflux','f4',('lat','lon',))

    latitudes[:] = lats
    longitudes[:] = lons

    latitudebs[:] = latbs
    longitudebs[:] = lonbs

    warmpool_array_netcdf[:] = warmpool_array

    warmpool_file.close()

