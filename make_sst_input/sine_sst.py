'''21/11/2017 Create a sinusoidally varying SST input based on Bordoni and Schneider 2010
   NB peaks will be at sin(lat) = sin(lat_max) -> Tmax = A + B (sin(23.5))^2
   Values at +- 90 will be A - B(1 +- 2sin(23.5))
'''

import numpy as np
import create_timeseries as cts
import xarray as xr
import matplotlib.pyplot as plt

# Set up details for output file
lons, lats, lonbs, latbs, nlon, nlat, nlonb, nlatb = cts.create_grid(manual_grid_option=False)
time_arr, day_number, ntime, time_units, time_bounds = cts.create_time_arr(num_years=1, is_climatology=True, time_spacing=360)

# Evaluate SST, getting max and min using formula above
timeangle = 23.5 * np.sin(np.pi/180. * np.arange(0.5,360.))

sst = np.zeros([360, 64])
for i in range(64):
    for j in range(360):
        sst[j,63-i] = 299. - 26. * (np.sin(lats[i]*np.pi/180.)**2 - 2.*np.sin(lats[i]*np.pi/180.) * np.sin (timeangle[j]*np.pi/180.))

print(np.min(np.min(sst)))
print(np.max(np.max(sst)))

#sst[sst<250.] = 250.

#plt.contourf(np.arange(0.5,360.), lats, sst.T, levels = np.arange(250.,305.,2.))
#plt.colorbar()
#plt.show()

# Put longitude back in!
data = np.repeat(np.expand_dims(sst, axis=2), 128, axis=2)

# Output to netcdf
number_dict={}
number_dict['nlat']=nlat
number_dict['nlon']=nlon
number_dict['nlatb']=nlatb
number_dict['nlonb']=nlonb
number_dict['ntime']=ntime

variable_name = 'sine_sst'
file_name = './sine_sst.nc'

cts.output_to_file(data, lats, lons, latbs, lonbs, None, None, time_arr, time_units, file_name, variable_name, number_dict)


for i in range(90,181,10):
    
    data_ss = data[i-1,:,:]

    data_ss = np.repeat(np.expand_dims(data_ss, axis=0), 360, axis=0)
    
    variable_name = 'sine_sst_' + str(i)
    file_name = './sine_sst_' + str(i) + '.nc'
    
    cts.output_to_file(data_ss, lats, lons, latbs, lonbs, None, None, time_arr, time_units, file_name, variable_name, number_dict)
    


