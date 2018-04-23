# -*- coding: utf-8 -*-s
import numpy as np
import create_timeseries as cts
import xarray as xr
import matplotlib.pyplot as plt


#lons, lats, lonbs, latbs, nlon, nlat, nlonb, nlatb = cts.create_grid(manual_grid_option=False)
time_arr, day_number, ntime, time_units, time_bounds = cts.create_time_arr(num_years=1, is_climatology=True, time_spacing=12)

data_file = 'sst_clim_amip.nc'
data = xr.open_dataset( data_file, decode_times=False)

lons = data.lon.values
lats = data.lat.values
lonbs = data.lonb.values
latbs = data.latb.values
nlon = len(lons)
nlat = len(lats)
nlonb = len(lonbs)
nlatb = len(latbs)

data = data.sst_clim_amip

# Take the zonal mean of the data
data_mean = data.mean('lon')
data_sym = data_mean
# Make sure DJF is the opposite of JJA etc. 
#data_sym = np.zeros(data_mean.values.shape)
#for i in range(0,6):
#    data_sym[i,:] = (data_mean[i,:].values + data_mean[i+6,::-1].values)/2.
#    data_sym[i+6,:] = data_sym[i,::-1]

# Put longitude back in!
data = np.repeat(np.expand_dims(data_sym, axis=2), 360, axis=2)

number_dict={}
number_dict['nlat']=nlat
number_dict['nlon']=nlon
number_dict['nlatb']=nlatb
number_dict['nlonb']=nlonb
number_dict['ntime']=ntime

variable_name = 'ap_amip_sst2'
file_name = './ap_amip_sst2.nc'

cts.output_to_file(data, lats, lons, latbs, lonbs, None, None, time_arr, time_units, file_name, variable_name, number_dict)




