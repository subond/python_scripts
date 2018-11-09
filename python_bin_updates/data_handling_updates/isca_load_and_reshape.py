'''13/06/2018 Load up '''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


def isca_load_and_reshape(run, field, months=[121,181], period_fac=1., filename='plev_pentad', days=False):
    '''Load and reshape multi-year data into year/day of year shape'''
    name_temp = '/disca/share/rg419/Data_moist/' + run + '/run%04d/' + filename + '.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    
    #read data into xarray 
    data = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    
    data = data[field]
    
    # Reshape to have dimensions ('year_no', 'xofyear', 'lat', 'lon')
    if days:
        data.coords['xofyear'] = np.mod( data.time, 360.*period_fac) + 0.5 
        year_no = np.repeat(np.arange(40.), int(360*period_fac))
        year_no = (data.time * 0.) + year_no[0: data.time.values.size]
        data.coords['year_no'] = year_no
        
    else:
        data.coords['xofyear'] = np.mod( data.time - 1., 360.*period_fac) //5 + 1.  
        year_no = np.repeat(np.arange(40.), int(72*period_fac))
        year_no = (data.time * 0.) + year_no[0: data.time.values.size]
        data.coords['year_no'] = year_no
    
    data_rs = data.set_index(time = ['year_no', 'xofyear'])
    data_rs = data_rs.unstack('time')
        
    if 'pfull' in data_rs.dims:
        data_rs = data_rs.transpose('year_no', 'xofyear', 'pfull', 'lat', 'lon')
    else:
        data_rs = data_rs.transpose('year_no', 'xofyear', 'lat', 'lon')
    
    data.close()
    data_rs.close()
    return data_rs
    

