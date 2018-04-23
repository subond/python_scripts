#check momentum budget closure of model
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from physics import mombudg_2d_an_fn
from data_handling import cell_area


for i in range(7,10):
    rundata = xr.open_dataset( '/scratch/rg419/Data_moist/co2_test_o3/run'+str(i)+'/atmos_monthly.nc',
      decode_times=False,  # no calendar so tell netcdf lib
      # choose how data will be broken down into manageable chunks.
       chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    u = rundata.ucomp.mean(('time','lon'))
    u.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(-10,60,5))
    plt.savefig('/scratch/rg419/plots/co2_climates/o3/u'+str(i)+'.png')
    plt.clf()
    
    t_surf = rundata.t_surf.mean(('time'))
    t_surf.plot.contourf(x='lon', y='lat', yincrease=False, levels = range(240,310,10))
    plt.savefig('/scratch/rg419/plots/co2_climates/o3/ts'+str(i)+'.png')
    plt.clf()    
    
    q = rundata.sphum.mean(('time','lon'))
    q.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,0.022,0.002))
    plt.savefig('/scratch/rg419/plots/co2_climates/o3/q'+str(i)+'.png')
    plt.clf()
    
    v = rundata.vcomp.mean(('time','lon'))
    v.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(-4,4.5,0.5))
    plt.savefig('/scratch/rg419/plots/co2_climates/o3/v'+str(i)+'.png')
    plt.clf()
    
    t = rundata.temp.mean(('time','lon'))
    t.plot.contourf(x='lat', y='pfull', yincrease=False, levels = range(180,305,10))
    plt.savefig('/scratch/rg419/plots/co2_climates/o3/t'+str(i)+'.png')
    plt.clf()
    
    
    tdt_rad = rundata.tdt_rad.mean(('time','lon'))*86400.
    tdt_rad.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(-5,5.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/o3/tdt_rad'+str(i)+'.png')
    plt.clf()
    
    tdt_solar = rundata.tdt_solar.mean(('time','lon'))*86400.
    tdt_solar.plot.contourf(x='lat', y='pfull', yincrease=False, levels = np.arange(0,2.1,0.25))
    plt.savefig('/scratch/rg419/plots/co2_climates/o3/tdt_solar'+str(i)+'.png')
    plt.clf()

