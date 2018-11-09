# 22/10/2018 Calculate the area mean of the wind speed to use as an input for a prescribed windspeed run

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh


def calc_wspeed_mean(run):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    coslat = np.cos(data.lat * np.pi/180.)
    v_area_weighted = data.wind_speed.mean('lon') * coslat # Weight u by coslat, to do area weighting
    v_mean = (v_area_weighted.sum('lat') / coslat.sum('lat')).mean('xofyear')
    
    print(v_mean)
    

if __name__ == "__main__":
    
    calc_wspeed_mean('sn_1.000_evap_fluxes')
    calc_wspeed_mean('sn_1.000_evap_fluxes_qflux')