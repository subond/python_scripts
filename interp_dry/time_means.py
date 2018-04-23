"""Script to take model output and produce a climatology. 
   23/11/2017 modified to allow files in hours to be converted to days using day_fac
   NB It's good to check that rundata.xofyear contains the values you'd expect"""

import numpy as np
import xarray as xr
import scipy.interpolate as spint

def time_means(run_fol, months, filename='atmos_daily', timeav='day', period_fac=1., write_netcdf=False, day_fac=1., offline_interp=None):
    #load data and average equivalent time periods over years
    if np.mod(months[0], 12*period_fac) != 1 and np.mod(months[1], 12*period_fac) != 1.:
        print 'WARNING: Non Jan/Dec start/end points'
    elif np.mod(months[0], 12*period_fac) != 1:
        print 'WARNING: Non Jan start point'
    elif np.mod(months[1], 12*period_fac) != 1:
        print 'WARNING: Non Dec end point'
    
    try:
        name_temp = '/scratch/rg419/Data_moist/' + run_fol + '/run%03d/'+filename+'.nc'
        names = [name_temp % m for m in range( months[0], months[1])  ]
    
        #print names
        #read data into xarray 
        rundata = xr.open_mfdataset( names,
            decode_times=False,  # no calendar so tell netcdf lib
            # choose how data will be broken down into manageable chunks.
            chunks={'time': 30})
    except:
        name_temp = '/scratch/rg419/Data_moist/' + run_fol + '/run%04d/'+filename+'.nc'
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
    if timeav == 'day':
        rundata.coords['xofyear'] = np.mod( rundata.time/day_fac, 360.*period_fac) + 0.5 
        data = rundata.groupby('xofyear').mean(('time'))
    elif timeav =='pentad':
        rundata.coords['xofyear'] = np.mod( rundata.time/day_fac - 1., 360.*period_fac) //5 + 1.  
        data = rundata.groupby('xofyear').mean(('time'))        
    elif timeav =='month':
        rundata.coords['xofyear'] = np.mod( rundata.time/day_fac, 360.*period_fac) //30 + 1        
        data = rundata.groupby('xofyear').mean(('time'))
    elif timeav =='season':
        rundata.coords['xofyear'] = np.mod( np.floor( (np.mod( rundata.time/day_fac, 360.*period_fac) //(30*period_fac) + 1.)/3. ) , 4);
        data = rundata.groupby('xofyear').mean(('time'))
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
            print data_var
            f = spint.interp1d(data.pfull, data[data_var], axis=1, fill_value='extrapolate', kind='quadratic')
            data_interp[data_var] = xr.DataArray(f(pfull_new), coords=[data.xofyear, pfull_new, data.lat, data.lon], dims=['xofyear', 'pfull', 'lat', 'lon'])
        
        data = data_interp
        
    if write_netcdf:
        filename = '/scratch/rg419/Data_moist/climatologies/'+run_fol+'.nc'
        data.to_netcdf(filename)
    
    print 'Data loaded'
    
    return data
    

if __name__ == "__main__":
    
    #offline_interp = ['div', 'height', 'omega', 'temp', 'teq', 'ucomp', 'vcomp', 'vor']
    #test=time_means('dry_zs', [9,21], filename='daily', timeav='pentad', write_netcdf=True, day_fac = 24., offline_interp=offline_interp)
    #test=time_means('dry_ep', [9,21], filename='daily', timeav='pentad', write_netcdf=True, day_fac = 24., offline_interp=offline_interp)
	
    test=time_means('sine_sst_10m', [121,349], filename='plev_pentad', timeav='pentad', write_netcdf=True)
	#test=time_means('ap10_co2', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)

	#test=time_means('1cont', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)
	#test=time_means('2cont', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)
	#test=time_means('sn_0.500', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True, period_fac=0.5)
	#test=time_means('sn_2.000', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True, period_fac=2.)
