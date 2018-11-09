'''28/09/2018 Combine NCEP-NCAR u-850 into one file and save'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
import pandas as pd
from climatology import bob_onset, scsm_onset, ism_onset


name_temp = '/disca/share/rg419/NCEP_NCAR/uwnd.%04d.nc'
names = [name_temp % m for m in range( 1948, 2017) ]
#read data into xarray 
data = xr.open_mfdataset(names)

#select 850 hPa level
data = data.sel(level=850.)

# output to netcdf on disca
data.to_netcdf('/disca/share/rg419/ncep_ncar_uwnd_daily_850_1948_2016.nc')
