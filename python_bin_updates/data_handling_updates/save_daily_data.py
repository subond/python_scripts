import numpy as np
import xarray as xr

def save_daily_data(run, months, var, filename='plev_daily', period_fac=1., lev=None):
    #load data and average equivalent time periods over years
    if np.mod(months[0], 12*period_fac) != 1 and np.mod(months[1], 12*period_fac) != 1.:
        print 'WARNING: Non Jan/Dec start/end points'
    elif np.mod(months[0], 12*period_fac) != 1:
        print 'WARNING: Non Jan start point'
    elif np.mod(months[1], 12*period_fac) != 1:
        print 'WARNING: Non Dec end point'
        
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
        
    #read data into xarray 
    rundata = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    
    print 'files opened, starting saving'
    
    if lev is not None:
        savevar = rundata[var].sel(pfull=lev).to_dataset()
        filename = '/scratch/rg419/Data_moist/climatologies/'+run+'_'+var+str(lev)+'_daily.nc'
    else:
        savevar = rundata[var].to_dataset()
        filename = '/scratch/rg419/Data_moist/climatologies/'+run+'_'+'_daily.nc'
    
    savevar = savevar.assign_coords(lonb=rundata.lonb)
    savevar = savevar.assign_coords(latb=rundata.latb)
    
    savevar.to_netcdf(filename, format='NETCDF3_CLASSIC')
    print 'data saved'
    


if __name__ == "__main__":


    #test = save_daily_data('qflux_5_300', [121,241], 'omega', lev=850.)
    #test = save_daily_data('qflux_10_300', [121,241], 'omega', lev=850.)
    #test = save_daily_data('qflux_15_300', [121,241], 'omega', lev=850.)
    #test = save_daily_data('qflux_20_300', [121,241], 'omega', lev=850.)
    #test = save_daily_data('qflux_25_300', [121,241], 'omega', lev=850.)
    #test = save_daily_data('ss_90.000', [121,181], 'omega', lev=850.)
    test = save_daily_data('ss_91.000', [121,181], 'omega', lev=850.)
    test = save_daily_data('ss_92.000', [121,181], 'omega', lev=850.)
    test = save_daily_data('ss_93.000', [121,181], 'omega', lev=850.)
    test = save_daily_data('ss_94.000', [121,181], 'omega', lev=850.)
    #test = save_daily_data('ss_95.000', [121,181], 'omega', lev=850.)
    #test = save_daily_data('ss_100.000', [121,181], 'omega', lev=850.)
    #test = save_daily_data('ss_105.000', [121,181], 'omega', lev=850.)
    

