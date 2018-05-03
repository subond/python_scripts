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
from data_handling_updates import flip_field


def p_cent_rate_max(run, days=None):
    # Get the maximum rate of change of the precipitation centroid latitude, and the latitude at which this occurs.
    
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%04d/plev_pentad.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    
    #read data into xarray 
    data = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    
    data = precip_load_and_reshape(data)
    
    # Locate precipitation centroid
    precip_centroid(data)
            
    # Get rate of movement of precip centroid
    if days[j]:
        dpcentdt = gr.ddt(data.p_cent, secperunit = 86400.) * 86400.
    else:
        dpcentdt = gr.ddt(data.p_cent) * 86400.
    
    pcent_max = data.p_cent.max('xofyear')            
    
    dpcentdt_ppos = dpcentdt.where(data.p_cent>=0.)   # Find precip centroid rate where precip centroid is in the northern hemisphere
    dpcentdt_max = dpcentdt_ppos.where(dpcentdt_ppos==dpcentdt_ppos.max(),drop=True)   # Find the maximum of the above
    if len(dpcentdt_max) > 1:
        dpcentdt_max = dpcentdt_max[0]
    pcent_dtmax = data.p_cent.sel(xofyear=dpcentdt_max.xofyear)    # Find the location of the preciptiation when the rate is maximum
        
    print(dpcentdt_max.values, pcent_dtmax.values)     # Print both
    
    max_rate.append(dpcentdt_max)
    max_rate_lat.append(pcent_dtmax)
    max_lat.append(pcent_max)

    return max_rate, max_rate_lat, max_lat



def precip_load_and_reshape(run, months=[121,181], period_fac=1.):
    
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
    
    #data_flip = flip_field(data_rs.precipitation)
    #data_flip = xr.Dataset({'precipitation': data_flip}, coords=data_flip.coords)
    
    # Get precipitation centroid
    #test = precip_centroid(data_rs)
    #test2= precip_centroid(data_flip)
    
    #for i in range(1):
    #    test.p_cent[i,].plot.line()
        #test2.p_cent[i,].plot.line()
    #plt.show()
    return data_rs
    
    
    
    # I want this to:
    # Load up all the years for a run
    # Rearrange them into year and pentad of year
    # Calculate the time and hemispheric mean precipitation centroid's max rate, lat of max rate, and max lat
    # Calculate the precipitation centroid's max rate, lat of max rate, and max lat for every year and both hemispheres
    # Firstly average this and compare with the value for the time and hemispheric mean
    # Then take the standard deviation of the above quantities (decide which mean to use at this point) 



if __name__ == "__main__":
    
    p_cent_error('sn_1.000')