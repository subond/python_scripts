#check momentum budget closure of model
from netCDF4 import Dataset
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from physics import mombudg_2d_an_fn
from data_handling import cell_area

rundata = xr.open_dataset( '/scratch/rg419/Data_moist/diag_tests/run1/atmos_daily.nc',
  decode_times=False,  # no calendar so tell netcdf lib
  # choose how data will be broken down into manageable chunks.
       chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

plt.figure(0)
uq = rundata.sphum_u.mean(('time','lon'))
uq.plot.contourf(x='lat', y='pfull', yincrease=False)

plt.figure(1)
u = rundata.ucomp.mean(('time','lon'))
u.plot.contourf(x='lat', y='pfull', yincrease=False)

plt.figure(2)
q = rundata.sphum.mean(('time','lon'))
q.plot.contourf(x='lat', y='pfull', yincrease=False)

plt.show()