"""
Load detrended era data and save climatologies of required parts
"""
import numpy as np
import xarray as xr
from physics import gradients as gr
#import matplotlib.pyplot as plt

filename_vo = '/scratch/rg419/obs_and_reanalysis/sep_levs_vo/era_vo_150_dt.nc'
filename_u = '/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_150_dt.nc'
filename_v = '/scratch/rg419/obs_and_reanalysis/sep_levs_v/era_v_150_dt.nc'

# Load vorticity
data_vo = xr.open_dataset(filename_vo, decode_times=False)

omega = 7.2921150e-5
f = 2.* omega * np.sin(data_vo.lat * np.pi/180.)

# Choose lon range
#lons = [data_vo.lon[i] for i in range(len(data_vo.lon)) if data_vo.lon[i] >= 60. and data_vo.lon[i] < 150.]

vor = data_vo['vo_150'] + f

# Take x and y gradients of vorticity
dvodx = gr.ddx(vor)
dvody = gr.ddy(vor, vector=False)
print 'gradients done'

# Load u and v
data_u = xr.open_dataset(filename_u, decode_times=False)
data_v = xr.open_dataset(filename_v, decode_times=False)


udvodx = (dvodx * data_u['u_150']).mean('year_no')
vdvody = (dvody * data_v['v_150']).mean('year_no')

#udvodx.sel(lon=lons).mean('lon').plot.contourf()
#plt.figure(2)
#vdvody.sel(lon=lons).mean('lon').plot.contourf()
#plt.figure(3)


print 'data ready, saving...'

ds = xr.Dataset({ 'udvodx' : (['lat', 'lon', 'day_of_yr'], udvodx)},
                 coords={'lat': ('lat', udvodx.lat),
                         'lon': ('lon', udvodx.lon),
                         'day_of_yr': ('day_of_yr', udvodx.day_of_yr)})
                         
filename = '/scratch/rg419/obs_and_reanalysis/era_udvodx_clim.nc'
ds.to_netcdf(filename)


ds = xr.Dataset({ 'vdvody' : (['lat', 'lon', 'day_of_yr'], vdvody)},
                 coords={'lat': ('lat', vdvody.lat),
                         'lon': ('lon', vdvody.lon),
                         'day_of_yr': ('day_of_yr', vdvody.day_of_yr)})
                         
filename = '/scratch/rg419/obs_and_reanalysis/era_vdvody_clim.nc'
ds.to_netcdf(filename)

ds = xr.Dataset({ 'horiz_adv' : (['lat', 'lon', 'day_of_yr'], (udvodx + vdvody))},
                 coords={'lat': ('lat', vdvody.lat),
                         'lon': ('lon', vdvody.lon),
                         'day_of_yr': ('day_of_yr', vdvody.day_of_yr)})
                         
filename = '/scratch/rg419/obs_and_reanalysis/era_horiz_adv_clim.nc'
ds.to_netcdf(filename)
