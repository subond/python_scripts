import numpy as np
import xarray as xr

def time_means(run_fol, months, filename='atmos_daily', timeav='day', period_fac=1., write_netcdf=False):
    #load data and average equivalent time periods over years
    if np.mod(months[0], 12*period_fac) != 1 and np.mod(months[1], 12*period_fac) != 1.:
        print 'WARNING: Non Jan/Dec start/end points'
    elif np.mod(months[0], 12*period_fac) != 1:
        print 'WARNING: Non Jan start point'
    elif np.mod(months[1], 12*period_fac) != 1:
        print 'WARNING: Non Dec end point'
    
    name_temp = '/scratch/rg419/Data_moist/' + run_fol + '/run%03d/'+filename+'.nc'
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
        rundata.coords['xofyear'] = np.mod( rundata.time, 360.*period_fac) + 0.5 
        data = rundata.groupby('xofyear').mean(('time'))
    elif timeav =='pentad':
        rundata.coords['xofyear'] = np.mod( rundata.time, 360.*period_fac) //5 + 1.  
        data = rundata.groupby('xofyear').mean(('time'))        
    elif timeav =='month':
        rundata.coords['xofyear'] = np.mod( rundata.time, 360.*period_fac) //30 + 1        
        data = rundata.groupby('xofyear').mean(('time'))
    elif timeav =='season':
        rundata.coords['xofyear'] = np.mod( np.floor( (np.mod( rundata.time, 360.*period_fac) //(30*period_fac) + 1.)/3. ) , 4);
        data = rundata.groupby('xofyear').mean(('time'))
    else:
        'invalid timeav'
        return
    
    if write_netcdf:
        filename = '/scratch/rg419/Data_moist/climatologies/'+run_fol+'.nc'
        data.to_netcdf(filename)
    
    print 'Data loaded'
    
    return data
    

if __name__ == "__main__":

	test=time_means('zs_sst', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)
	test=time_means('ap10_qflux', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)
	test=time_means('ap10_co2', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)

	#test=time_means('1cont', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)
	#test=time_means('2cont', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)
	#test=time_means('sn_0.500', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True, period_fac=0.5)
	#test=time_means('sn_2.000', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True, period_fac=2.)
