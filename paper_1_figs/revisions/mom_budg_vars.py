"""
Load variables for momentum budget and save in one dataset for ease of use
"""
import numpy as np
import xarray as xr

dir = '/scratch/rg419/obs_and_reanalysis/'

# Combine multi-level climatologies
lev_list = range(50,350,50)
no_levs = len(lev_list)

u_in = xr.open_dataset(dir + 'sep_levs_u/era_u_150_clim.nc', decode_times=False)
u_shape =  u_in['u_150'].values.shape

u = np.zeros((u_shape[2], no_levs, u_shape[0], u_shape[1]))
w = np.zeros((u_shape[2], no_levs, u_shape[0], u_shape[1]))
uw = np.zeros((u_shape[2], no_levs, u_shape[0], u_shape[1]))

i = 0
for lev in lev_list:
    u_in = xr.open_dataset(dir + 'sep_levs_u/era_u_' + str(lev) +'_clim.nc', decode_times=False)
    uw_in = xr.open_dataset(dir + 'era_uw_' + str(lev) +'_clim.nc', decode_times=False)
    w_in = xr.open_dataset(dir + 'sep_levs_w/era_w_' + str(lev) +'_clim.nc', decode_times=False)
    
    u[:,i,:,:] = u_in.transpose('day_of_yr', 'lat', 'lon')['u_' + str(lev)].values
    w[:,i,:,:] = w_in.transpose('day_of_yr', 'lat', 'lon')['w_' + str(lev)].values
    uw[:,i,:,:] = uw_in.transpose('day_of_yr', 'lat', 'lon')['uw_' + str(lev)].values
    
    i = i+1


# Load single level climatologies
v_clim = xr.open_dataset(dir + 'era_v_clim.nc', decode_times=False)
v_clim = v_clim.transpose('day_of_yr', 'lat', 'lon')['v']
uu_clim = xr.open_dataset(dir + 'era_uu_150_clim.nc', decode_times=False)
uu_clim = uu_clim.transpose('day_of_yr', 'lat', 'lon')['uu_150']
uv_clim = xr.open_dataset(dir + 'era_uv_150_clim.nc', decode_times=False)
uv_clim = uv_clim.transpose('day_of_yr', 'lat', 'lon')['uv_150']

z_clim = xr.open_dataset(dir + 'era_z_clim.nc', decode_times=False)
z_clim = z_clim.transpose('day_of_yr', 'lat', 'lon')['z']


# Write all to netcdf
filename = '/scratch/rg419/obs_and_reanalysis/era_mom_vars.nc'
ds = xr.Dataset({ 'ucomp': (['day_of_yr', 'pfull', 'lat', 'lon'], u),
                  'vcomp': (['day_of_yr', 'lat', 'lon'], v_clim),
                  'omega': (['day_of_yr', 'pfull', 'lat', 'lon'], w),
                  'ucomp_sq': (['day_of_yr', 'lat', 'lon'], uu_clim),
                  'ucomp_vcomp': (['day_of_yr', 'lat', 'lon'], uv_clim),
                  'ucomp_omega': (['day_of_yr', 'pfull', 'lat', 'lon'], uw),
                  'phi': (['day_of_yr', 'lat', 'lon'], z_clim)},
                    coords={'lat': ('lat', u_in.lat),
                         'lon': ('lon', u_in.lon),
                         'day_of_yr': ('day_of_yr', u_in.day_of_yr),
                          'pfull': ('pfull', np.arange(50.,350.,50.))})

ds.to_netcdf(filename)

