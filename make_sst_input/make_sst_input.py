# -*- coding: utf-8 -*-s
import numpy as np
import create_timeseries as cts
import xarray as xr
import matplotlib.pyplot as plt


lons, lats, lonbs, latbs, nlon, nlat, nlonb, nlatb = cts.create_grid(manual_grid_option=False)
time_arr, day_number, ntime, time_units, time_bounds = cts.create_time_arr(num_years=1, is_climatology=True, time_spacing=72)

data_file = '/scratch/rg419/Data_moist/climatologies/sn_1.000.nc'
data = xr.open_dataset( data_file, decode_times=False)
data = data.t_surf

# Take the zonal mean of the data
data_mean = data.mean('lon')

# Make sure DJF is the opposite of JJA etc. 
data_sym = np.zeros(data_mean.values.shape)
for i in range(0,36):
    data_sym[i,:] = (data_mean[i,:].values + data_mean[i+36,::-1].values)/2.
    data_sym[i+36,:] = data_sym[i,::-1]

# Put longitude back in!
data = np.repeat(np.expand_dims(data_sym, axis=2), 128, axis=2)


number_dict={}
number_dict['nlat']=nlat
number_dict['nlon']=nlon
number_dict['nlatb']=nlatb
number_dict['nlonb']=nlonb
number_dict['ntime']=ntime

variable_name = 'sn_1.000_sst'
file_name = './sn_1.000_sst.nc'

cts.output_to_file(data, lats, lons, latbs, lonbs, None, None, time_arr, time_units, file_name, variable_name, number_dict)


time_arr, day_number, ntime, time_units, time_bounds = cts.create_time_arr(num_years=1, is_climatology=True, time_spacing=360)

for i in range(36,55,2):
    
    data_ss = (data[i-1,:,:] + data[i,:,:])/2.

    data_ss = np.repeat(np.expand_dims(data_ss, axis=0), 360, axis=0)
    
    variable_name = 'sn_1.000_sst_' + str(i)
    file_name = './sn_1.000_sst_' + str(i) + '.nc'
    
    cts.output_to_file(data_ss, lats, lons, latbs, lonbs, None, None, time_arr, time_units, file_name, variable_name, number_dict)
    


