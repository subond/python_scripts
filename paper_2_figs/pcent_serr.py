"""
Evaluate standard error of precipitation centroid 1/05/2018
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_centroid
from data_handling_updates import gradients as gr
from pylab import rcParams
from data_handling_updates import make_sym


def p_cent_error(run, months=[121,181], period_fac=1.):
    
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%04d/plev_pentad.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    
    #read data into xarray 
    data = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    
    data = data.precipitation
    
    # Reshape to have dimensions ('year_no', 'xofyear', 'lat', 'lon')
    data.coords['xofyear'] = np.mod( data.time - 1., 360.*period_fac) //5 + 1.  
    
    year_no = np.repeat(np.arange(40.), 72)
    year_no = (data.time * 0.) + year_no[0: data.time.values.size]
    data.coords['year_no'] = year_no
    
    data_rs = data.set_index(time = ['year_no', 'xofyear'])
    data_rs = data_rs.unstack('time')
    
    data_rs = xr.Dataset({'precipitation': data_rs}, coords=data_rs.coords)
    data_rs = data_rs.transpose('year_no', 'xofyear', 'lat', 'lon')
    
    # Get precipitation centroid
    test = precip_centroid(data_rs)
    
    for i in range(5):
        test.p_cent[i,].plot.line()
    plt.show()



if __name__ == "__main__":
    
    p_cent_error('sn_1.000')