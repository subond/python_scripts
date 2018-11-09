"""
Plot angular momentum conserving temperature for different runs (06/09/2018)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, model_constants as mc
from climatology import precip_centroid
from pylab import rcParams
from hadley_cell import mass_streamfunction


def get_latmax(data_in, max_by=[]):
    if len(max_by) == 0.:
        data_max = data_in.max('lat')
        data_max_lat = data_in.lat.values[data_in.argmax('lat').values]
        data_max_lat = xr.DataArray(data_max_lat, coords=[data_in.xofyear], dims=['xofyear'])
    else:
        data_max_lat = max_by.lat.values[max_by.argmax('lat').values]
        data_max = np.zeros(len(data_in.xofyear.values))
        for i in range(len(data_in.xofyear.values)):
            data_max[i] = data_in[i,:].sel(lat=data_max_lat[i]).values
        data_max = xr.DataArray(data_max, coords=[data_in.xofyear], dims=['xofyear'])
        data_max_lat = xr.DataArray(data_max_lat, coords=[data_in.xofyear], dims=['xofyear'])
    return data_max, data_max_lat


def plot_ang_mom_temp(run, pentad, lev=850., color='k', use_omega=False, rotfac=1.):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    convTtotheta=(1000./data.pfull)**mc.kappa
    theta = data.temp*convTtotheta
    theta_lev = theta.sel(pfull=lev).mean('lon')
    
    omega_max, omega_max_lat = get_latmax(-1.*data.omega.sel(pfull=lev).mean('lon'))
    
    if use_omega:
        theta_max, theta_max_lat = get_latmax(theta_lev, max_by=-1.*data.omega.sel(pfull=lev).mean('lon'))
    else:
        theta_max, theta_max_lat = get_latmax(theta_lev)
    
    sinphi = np.sin(data.lat * np.pi/180.)
    cosphi = np.cos(data.lat * np.pi/180.)
    sinphimax = np.sin(theta_max_lat * np.pi/180.)
    
    theta_m = theta_max - 300.*(mc.omega*rotfac)**2.*mc.a**2./(2.*mc.grav*14000.) * (sinphi**2.-sinphimax**2.)**2./cosphi**2.
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -60. and data.lat[i] <= 60.]
    
    theta_lev.sel(xofyear=pentad, lat=lats).plot(color=color)
    theta_m.sel(xofyear=pentad, lat=lats).plot(color=color, linestyle='--')
    #plt.plot([omega_max_lat.sel(xofyear=pentad),omega_max_lat.sel(xofyear=pentad)],[295,310],color=color,linestyle=':')
    plt.ylim([280,320])

plot_ang_mom_temp('sn_1.000', 55, color='C0', use_omega=False)
plot_ang_mom_temp('rt_0.750', 55, color='C1', use_omega=False)
plot_ang_mom_temp('rt_1.250', 55, color='C1', use_omega=False)
plot_ang_mom_temp('mld_20', 55, color='C2', use_omega=False)
plot_ang_mom_temp('mld_2.5', 55, color='C2', use_omega=False)
plt.show()