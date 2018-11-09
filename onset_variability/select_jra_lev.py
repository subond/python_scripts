'''28/09/2018 Combine NCEP-NCAR u-850 into one file and save'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
import pandas as pd
from climatology import bob_onset, scsm_onset, ism_onset


data = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/ucomp_daily/atmos_daily_together.nc')
data = data.sel(lev=20000.)
data.to_netcdf('/disca/share/rg419/jra_ucomp_daily_200.nc')

data.close()

data = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_daily/atmos_daily_together.nc')
data = data.sel(lev=20000.)
data.to_netcdf('/disca/share/rg419/jra_vcomp_daily_200.nc')

data.close()