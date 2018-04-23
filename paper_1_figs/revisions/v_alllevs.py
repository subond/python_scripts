"""
Load variables for momentum budget and save in one dataset for ease of use
"""
import numpy as np
import xarray as xr

dir = '/scratch/rg419/obs_and_reanalysis/'

# Combine multi-level climatologies
lev_list = range(1000,0,-50)
no_levs = len(lev_list)

v_in = xr.open_dataset(dir + 'sep_levs_v/era_v_150_clim.nc', decode_times=False)
v_shape =  v_in['v_150'].values.shape

v = np.zeros((v_shape[2], no_levs, v_shape[0], v_shape[1]))


i = 0
for lev in lev_list:
    v_in = xr.open_dataset(dir + 'sep_levs_v/era_v_' + str(lev) +'_clim.nc', decode_times=False)
    
    v[:,i,:,:] = v_in.transpose('day_of_yr', 'lat', 'lon')['v_' + str(lev)].values
    
    i = i+1


# Write all to netcdf
filename = '/scratch/rg419/obs_and_reanalysis/era_v_alllevs.nc'
ds = xr.Dataset({ 'vcomp': (['day_of_yr', 'pfull', 'lat', 'lon'], v)},
                    coords={'lat': ('lat', v_in.lat),
                         'lon': ('lon', v_in.lon),
                         'day_of_yr': ('day_of_yr', v_in.day_of_yr),
                          'pfull': ('pfull', np.arange(1000.,0.,-50.))})

ds.to_netcdf(filename)

