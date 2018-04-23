#check momentum budget closure of model
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from physics import mombudg_2d_an_fn
from data_handling import cell_area


for i in range(49,61):
    rundata = xr.open_dataset( '/scratch/rg419/Data_moist/co2_test_frierson/run'+str(i)+'/atmos_monthly.nc',
      decode_times=False,  # no calendar so tell netcdf lib
      # choose how data will be broken down into manageable chunks.
       chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    u = rundata.ucomp.mean(('time','lon'))
    u.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(-10,60,5))
    plt.savefig('/scratch/rg419/plots/co2_climates/frierson/u'+str(i)+'.png')
    plt.clf()
    
    tdt_rad = rundata.tdt_rad.mean(('time','lon'))*86400.
    tdt_rad.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(-5,5.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/frierson/tdt_rad'+str(i)+'.png')
    plt.clf()
    
    tdt_solar = rundata.tdt_solar.mean(('time','lon'))*86400.
    tdt_solar.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,2.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/frierson/tdt_solar'+str(i)+'.png')
    plt.clf()
    
    t = rundata.temp.mean(('time','lon'))
    t.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(180,305,10))
    plt.savefig('/scratch/rg419/plots/co2_climates/frierson/t'+str(i)+'.png')
    plt.clf()
    
    q = rundata.sphum.mean(('time','lon'))
    q.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,0.036,0.002))
    plt.savefig('/scratch/rg419/plots/co2_climates/frierson/q'+str(i)+'.png')
    plt.clf()


    rundata = xr.open_dataset( '/scratch/rg419/Data_moist/co2_test_rrtm/run'+str(i)+'/atmos_monthly.nc',
      decode_times=False,  # no calendar so tell netcdf lib
      # choose how data will be broken down into manageable chunks.
       chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    u = rundata.ucomp.mean(('time','lon'))
    u.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(-10,60,5))
    plt.savefig('/scratch/rg419/plots/co2_climates/rrtm/u'+str(i)+'.png')
    plt.clf()
    
    tdt_rad = rundata.tdt_rad.mean(('time','lon'))*86400.
    tdt_rad.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(-5,5.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/rrtm/tdt_rad'+str(i)+'.png')
    plt.clf()
    
    tdt_solar = rundata.tdt_sw.mean(('time','lon'))*86400.
    tdt_solar.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,2.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/rrtm/tdt_solar'+str(i)+'.png')
    plt.clf()
    
    t = rundata.temp.mean(('time','lon'))
    t.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(180,305,10))
    plt.savefig('/scratch/rg419/plots/co2_climates/rrtm/t'+str(i)+'.png')
    plt.clf()
    
    q = rundata.sphum.mean(('time','lon'))
    q.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,0.036,0.002))
    plt.savefig('/scratch/rg419/plots/co2_climates/rrtm/q'+str(i)+'.png')
    plt.clf()



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
    
    
for i in range(25,37):

    rundata = xr.open_dataset( '/scratch/rg419/Data_moist/co2_test_1/run'+str(i)+'/atmos_monthly.nc',
      decode_times=False,  # no calendar so tell netcdf lib
      # choose how data will be broken down into manageable chunks.
       chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    u = rundata.ucomp.mean(('time','lon'))
    u.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(-10,60,5))
    plt.savefig('/scratch/rg419/plots/co2_climates/me/u'+str(i)+'.png')
    plt.clf()
    
    tdt_rad = rundata.tdt_rad.mean(('time','lon'))*86400.
    tdt_rad.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(-5,5.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/me/tdt_rad'+str(i)+'.png')
    plt.clf()
    
    tdt_solar = rundata.tdt_solar.mean(('time','lon'))*86400.
    tdt_solar.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,2.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/me/tdt_solar'+str(i)+'.png')
    plt.clf()
    
    t = rundata.temp.mean(('time','lon'))
    t.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(180,305,10))
    plt.savefig('/scratch/rg419/plots/co2_climates/me/t'+str(i)+'.png')
    plt.clf()
    
    q = rundata.sphum.mean(('time','lon'))
    q.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,0.036,0.002))
    plt.savefig('/scratch/rg419/plots/co2_climates/me/q'+str(i)+'.png')
    plt.clf()