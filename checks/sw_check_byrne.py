#check momentum budget closure of model
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from physics import mombudg_2d_an_fn
from data_handling import cell_area


for i in range(1,6):

    rundata = xr.open_dataset( '/scratch/rg419/Data_moist/co2_test_byrne/run'+str(i)+'/atmos_monthly.nc',
      decode_times=False,  # no calendar so tell netcdf lib
      # choose how data will be broken down into manageable chunks.
       chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    tdt_solar = rundata.tdt_solar.mean(('time','lon'))*86400.
    tdt_solar.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,2.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrne/tdt_solar'+str(i)+'.png')
    plt.clf()

    
    