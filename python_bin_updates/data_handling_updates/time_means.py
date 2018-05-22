"""Script to take model output and produce a climatology. 
   23/11/2017 modified to allow files in hours to be converted to days using day_fac
   NB It's good to check that rundata.xofyear contains the values you'd expect"""

import numpy as np
import xarray as xr
import scipy.interpolate as spint

def time_means(run_fol, months, filename='atmos_daily', timeav='day', period_fac=1., write_netcdf=False, day_fac=1., offline_interp=None, name_out=None):
    #load data and average equivalent time periods over years
    if np.mod(months[0], 12*period_fac) != 1 and np.mod(months[1], 12*period_fac) != 1.:
        print ('WARNING: Non Jan/Dec start/end points')
    elif np.mod(months[0], 12*period_fac) != 1:
        print ('WARNING: Non Jan start point')
    elif np.mod(months[1], 12*period_fac) != 1:
        print ('WARNING: Non Dec end point')
    
    try:
        name_temp = '/disca/share/rg419/Data_moist/' + run_fol + '/run%03d/'+filename+'.nc'
        names = [name_temp % m for m in range( months[0], months[1])  ]
    
        #print names
        #read data into xarray 
        rundata = xr.open_mfdataset( names,
            decode_times=False,  # no calendar so tell netcdf lib
            # choose how data will be broken down into manageable chunks.
            chunks={'time': 30})
    except:
        name_temp = '/disca/share/rg419/Data_moist/' + run_fol + '/run%04d/'+filename+'.nc'
        names = [name_temp % m for m in range( months[0], months[1])  ]
    
        #print names
        #read data into xarray 
        rundata = xr.open_mfdataset( names,
            decode_times=False,  # no calendar so tell netcdf lib
            # choose how data will be broken down into manageable chunks.
            chunks={'time': 30})


    if timeav == '6hour':
        rundata.coords['xofyear'] = np.mod( rundata.time, 8640.*period_fac) 
        data = rundata.groupby('xofyear').mean(('time'))   
        filename = '/disca/share/rg419/Data_moist/climatologies/'+run_fol+'_6hour.nc'
    if timeav == 'day':
        #rundata.coords['xofyear'] = np.mod( rundata.time/day_fac -1., 360.*period_fac) + 0.5 
        rundata.coords['xofyear'] = np.mod( rundata.time/day_fac , 360.*period_fac) + 0.5 
        print (rundata.xofyear.values)
        data = rundata.groupby('xofyear').mean(('time'))
        filename = '/disca/share/rg419/Data_moist/climatologies/'+run_fol+'_day.nc'
    elif timeav =='pentad':
        rundata.coords['xofyear'] = np.mod( rundata.time/day_fac - 1., 360.*period_fac) //5 + 1.  
        print (rundata.xofyear.values)
        data = rundata.groupby('xofyear').mean(('time'))        
        filename = '/disca/share/rg419/Data_moist/climatologies/'+run_fol+'.nc'
    elif timeav =='month':
        rundata.coords['xofyear'] = np.mod( rundata.time/day_fac, 360.*period_fac) //30 + 1     
        #print rundata.xofyear.values   
        data = rundata.groupby('xofyear').mean(('time'))
        filename = '/disca/share/rg419/Data_moist/climatologies/'+run_fol+'_month.nc'
    elif timeav =='season':
        rundata.coords['xofyear'] = np.mod( np.floor( (np.mod( rundata.time/day_fac, 360.*period_fac) //(30*period_fac) + 1.)/3. ) , 4);
        data = rundata.groupby('xofyear').mean(('time'))
        filename = '/disca/share/rg419/Data_moist/climatologies/'+run_fol+'_season.nc'
    else:
        'invalid timeav'
        return
    
    
    if not offline_interp == None:
        pfull_new = np.arange(1000., 0., -50.)        
        data_interp = xr.Dataset({'ps': (['xofyear', 'lat', 'lon'], data.ps),
                         'h_trop':  (['xofyear', 'lat', 'lon'], data.h_trop)},
                         coords={'xofyear': ('xofyear', data.xofyear),
                                 'pfull': ('pfull', pfull_new),
                                   'lat': ('lat', data.lat),
                                   'lon': ('lon', data.lon),
                                   'latb': ('latb', data.latb),
                                   'lonb': ('lonb', data.lonb)})
        for data_var in offline_interp:
            print (data_var)
            f = spint.interp1d(data.pfull, data[data_var], axis=1, fill_value='extrapolate', kind='quadratic')
            data_interp[data_var] = xr.DataArray(f(pfull_new), coords=[data.xofyear, pfull_new, data.lat, data.lon], dims=['xofyear', 'pfull', 'lat', 'lon'])
        
        data = data_interp
        
    if write_netcdf:
        #filename = '/disca/share/rg419/Data_moist/climatologies/'+run_fol+'_day.nc'
        if name_out == None:
            data.to_netcdf(filename)
        else:
            data.to_netcdf(name_out)
    
    print ('Data loaded')
    
    return data
    

if __name__ == "__main__":
    
    
    
    test=time_means('sn_0.250', [121,211], filename='plev_daily_mean', timeav='day', write_netcdf=True, period_fac=0.25)


