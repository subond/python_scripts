#check momentum budget closure of model
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from physics import mombudg_2d_an_fn
from data_handling import cell_area


for i in range(60,61):

    rundata = xr.open_dataset( '/scratch/rg419/Data_moist/co2_test_byrne/run'+str(i)+'/atmos_monthly.nc',
      decode_times=False,  # no calendar so tell netcdf lib
      # choose how data will be broken down into manageable chunks.
       chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    u = rundata.ucomp.mean(('time','lon'))
    u.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(-10,60,5))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrne/u'+str(i)+'.png')
    plt.clf()
    
    tdt_rad = rundata.tdt_rad.mean(('time','lon'))*86400.
    tdt_rad.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(-5,5.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrne/tdt_rad'+str(i)+'.png')
    plt.clf()
    
    tdt_solar = rundata.tdt_solar.mean(('time','lon'))*86400.
    tdt_solar.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,2.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrne/tdt_solar'+str(i)+'.png')
    plt.clf()
    
    t = rundata.temp.mean(('time','lon'))
    t.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(180,305,10))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrne/t'+str(i)+'.png')
    plt.clf()
    
    q = rundata.sphum.mean(('time','lon'))
    q.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,0.036,0.002))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrne/q'+str(i)+'.png')
    plt.clf()
    
    
    rundata = xr.open_dataset( '/scratch/rg419/Data_moist/co2_test_byrnemod/run'+str(i)+'/atmos_monthly.nc',
      decode_times=False,  # no calendar so tell netcdf lib
      # choose how data will be broken down into manageable chunks.
       chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    ub = rundata.ucomp.mean(('time','lon'))
    ub.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(-10,60,5))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrnemod/u'+str(i)+'.png')
    plt.clf()
    
    tdt_rad = rundata.tdt_rad.mean(('time','lon'))*86400.
    tdt_rad.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(-5,5.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrnemod/tdt_rad'+str(i)+'.png')
    plt.clf()
    
    tdt_solar = rundata.tdt_solar.mean(('time','lon'))*86400.
    tdt_solar.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,2.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrnemod/tdt_solar'+str(i)+'.png')
    plt.clf()
    
    tb = rundata.temp.mean(('time','lon'))
    tb.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(180,305,10))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrnemod/t'+str(i)+'.png')
    plt.clf()
    
    q = rundata.sphum.mean(('time','lon'))
    q.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,0.036,0.002))
    plt.savefig('/scratch/rg419/plots/co2_climates/byrnemod/q'+str(i)+'.png')
    plt.clf()
    
    (ub-u).plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(-10,11,1))
    u.plot.contour(x='lat', y='pfull', yincrease=False, colors='k', levels = range(-10,60,5), add_colorbar=False)
    
    plt.figure(2)
    (tb-t).plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(-10,11,1))
    t.plot.contour(x='lat', y='pfull', yincrease=False, colors='k', levels = range(180,305,10), add_colorbar=False)
    
    plt.show()
    