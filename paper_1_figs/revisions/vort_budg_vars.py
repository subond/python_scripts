"""
Load variables for vorticity budget and save in one dataset for ease of use
"""
import numpy as np
import xarray as xr
from physics import gradients as gr
import matplotlib.pyplot as plt

dir = '/scratch/rg419/obs_and_reanalysis/'

# Load single level climatologies
u_clim = xr.open_dataset(dir + 'sep_levs_u/era_u_150_clim.nc', decode_times=False)
u_clim = u_clim.transpose('day_of_yr', 'lat', 'lon')['u_150']

v_clim = xr.open_dataset(dir + 'sep_levs_v/era_v_150_clim.nc', decode_times=False)
v_clim = v_clim.transpose('day_of_yr', 'lat', 'lon')['v_150']

vo_clim = xr.open_dataset(dir + 'sep_levs_vo/era_vo_150_clim.nc', decode_times=False)
vo_clim = vo_clim.transpose('day_of_yr', 'lat', 'lon')['vo_150']

div_clim = xr.open_dataset(dir + 'processed_era/era_d_clim.nc', decode_times=False)
div_clim = div_clim.transpose('day_of_yr', 'lat', 'lon')['d']

stretching_clim = xr.open_dataset(dir + 'processed_era/era_vod_150_clim.nc', decode_times=False)
stretching_clim = stretching_clim.transpose('day_of_yr', 'lat', 'lon')['vod_150']

horiz_adv_clim = xr.open_dataset(dir + 'processed_era/era_horiz_adv_clim.nc', decode_times=False)
horiz_adv_clim = horiz_adv_clim.transpose('day_of_yr', 'lat', 'lon')['horiz_adv']

omega = 7.2921150e-5
f = 2.* omega * np.sin(v_clim.lat * np.pi/180.)

stretching_clim = -1. * (stretching_clim + div_clim * f)

horiz_adv_clim = -1. * horiz_adv_clim

# Write all to netcdf
filename = '/scratch/rg419/obs_and_reanalysis/era_vort_vars.nc'
ds = xr.Dataset({ 'ucomp': (['day_of_yr', 'lat', 'lon'], u_clim),
                  'vcomp': (['day_of_yr', 'lat', 'lon'], v_clim),
                  'vor': (['day_of_yr', 'lat', 'lon'], vo_clim + f),
                  'div': (['day_of_yr', 'lat', 'lon'], div_clim),
                  'horiz_adv': (['day_of_yr', 'lat', 'lon'], horiz_adv_clim),
                  'stretching': (['day_of_yr', 'lat', 'lon'], stretching_clim)},
                    coords={'lat': ('lat', u_clim.lat),
                         'lon': ('lon', u_clim.lon),
                         'day_of_yr': ('day_of_yr', u_clim.day_of_yr)})

ds.to_netcdf(filename)

