'''Open files for last 2 years of spin-up and last 2 years of perturbed equilbration test runs
   Make SSTs to test the atmosphere response when the change is immediate'''

# -*- coding: utf-8 -*-s
import numpy as np
import create_timeseries as cts
import xarray as xr
import matplotlib.pyplot as plt


lons, lats, lonbs, latbs, nlon, nlat, nlonb, nlatb = cts.create_grid(manual_grid_option=False)
time_arr, day_number, ntime, time_units, time_bounds = cts.create_time_arr(num_years=1, is_climatology=True, time_spacing=360)
        
name_temp = '/scratch/rg419/Data_moist/ss_eq_2.5/run%04d/plev_monthly.nc'
names = [name_temp % m for m in range(49,61)  ]
data_eq = xr.open_mfdataset( names, decode_times=False, chunks={'time': 30})
data_eq = data_eq.t_surf

name_temp = '/scratch/rg419/Data_moist/ss_eq_2.5/run%04d/plev_pentad.nc'
names = [name_temp % m for m in range(109,121)  ]
data_sol = xr.open_mfdataset( names, decode_times=False, chunks={'time': 30})
data_sol = data_sol.t_surf


# Take the zonal mean of the data
data_eq = data_eq.mean(('lon','time'))
data_sol = data_sol.mean(('lon','time'))

# Take hemispheric mean for data_eq
data_eq = (data_eq[::] + data_eq[::-1])/2.

check = np.expand_dims(data_eq, axis=2)

# Put time and longitude back in
data_eq = np.repeat(np.expand_dims(data_eq, axis=0), 360, axis=0)
data_eq = np.repeat(np.expand_dims(data_eq, axis=2), 128, axis=2)
data_sol = np.repeat(np.expand_dims(data_sol, axis=0), 360, axis=0)
data_sol = np.repeat(np.expand_dims(data_sol, axis=2), 128, axis=2)

number_dict={}
number_dict['nlat']=nlat
number_dict['nlon']=nlon
number_dict['nlatb']=nlatb
number_dict['nlonb']=nlonb
number_dict['ntime']=ntime

variable_name = 'eq_sst'
file_name = './eq_sst.nc'
cts.output_to_file(data_eq, lats, lons, latbs, lonbs, None, None, time_arr, time_units, file_name, variable_name, number_dict)

variable_name = 'sol_sst'
file_name = './sol_sst.nc'
cts.output_to_file(data_sol, lats, lons, latbs, lonbs, None, None, time_arr, time_units, file_name, variable_name, number_dict)


