#Create ozone single time file.
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

a= 6376.0e3 #radius used in model

name='/scratch/rg419/GFDL_model/GFDLmoistModel/exp/q_flux_test/input/ocean_qflux.nc'
rundata = xr.open_dataset( name,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

#rundata.ocean_qflux.values[rundata.ocean_qflux.values==0] = np.nan

def integrate_lat(data):
    rundata.coords['dlat'] = ('lat',np.diff(data.latb))
    rundata.coords['dlon'] = ('lon',np.diff(data.lonb))
    integrand = rundata.ocean_qflux * rundata.dlat*np.pi/180. * a * a * np.cos(rundata.lat * np.pi/180.) * rundata.dlon*np.pi/180.
    
    integrand_test = (integrand).sum(('lon'))

    rundata['heat_trans'] = (('time','lat','lon'), np.cumsum(integrand.values, axis=1) )

    rundata['heat_trans_test'] = (('time','lat'), np.cumsum(integrand_test.values, axis=1) )
    
    rundata.heat_trans.values[rundata.ocean_qflux.values==0] = np.nan
    

    #rundata.heat_trans_test.mean('time').plot()
    rundata['heat_trans_mean'] = (('time','lat'),np.nanmean(rundata.heat_trans.values, axis=2)) #.nanmean(('time','lon'))#.plot()
    rundata.heat_trans_mean[0,:].plot()
    
    #rundata.heat_trans[0,:,:].mean(('lon')).plot()
    #rundata.heat_trans[0,:,20].plot()

    #plt.figure(2)
    #rundata.heat_trans[0,63,:].plot()

    plt.show()
    
integrate_lat(rundata)

