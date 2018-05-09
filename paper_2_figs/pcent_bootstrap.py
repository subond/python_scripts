"""
Use a bootstrapping method to estimate the distribution of the precipitation centroid's max rate, lat of max rate, and max lat 3/5/2018
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_centroid
from data_handling_updates import flip_field, gradients as gr
from pylab import rcParams
from pcent_rate_max import p_cent_rate_max, precip_load_and_reshape
from random import randint


    # I want this to:
    # Load up all the years for a run
    # Rearrange them into year and pentad of year
    # Calculate the time and hemispheric mean precipitation centroid's max rate, lat of max rate, and max lat
    # Calculate the precipitation centroid's max rate, lat of max rate, and max lat for every year and both hemispheres
    # Firstly average this and compare with the value for the time and hemispheric mean
    # Then take the standard deviation of the above quantities (decide which mean to use at this point) 

def bootstrap_sample(size):
    # Get the sample for the bootstrap
    sample=[]
    for i in range(size):
        a = randint(0,size-1)
        sample.append(a)
    return sample


def get_data(run, months=[121,481], period_fac=1.):
    # Load in the data, reshape it, flip it and add flipped data onto the end
    # This way, all means will be hemispheric and time, and both hemispheres will be included in sampling
    print(run)
    data = precip_load_and_reshape(run, months=months, period_fac=period_fac)
    
    n_years = float(data.year_no.size)
    
    precip_sh = flip_field(data.precipitation)
    precip_sh = xr.DataArray(precip_sh.values, coords=[n_years+np.arange(n_years), data.xofyear, data.lat, data.lon], dims=['year_no','xofyear','lat','lon'])
    
    data_sh = xr.Dataset({'precipitation': precip_sh}, coords=precip_sh.coords)
    data_extended = xr.concat([data, data_sh], 'year_no')
    
    
    return data_extended


def bootstrap_method(data, days=False):
    
    max_rate = []
    max_rate_lat = []
    max_lat = []
    
    for boot in range(1000):
        print(boot)
    
        sample = bootstrap_sample(data.year_no.size)
    
        precip_sample = data.sel(year_no=sample).mean('year_no')
        precip_centroid(precip_sample)
    
        # Next step is to unpick pcent_rate_max function so that you can just hand it the data for the sample mean
        # Get the numbers for 1000 or so samples, find 5 and 95 percentile  
     
        if days:
            dpcentdt = gr.ddt(precip_sample.p_cent, secperunit = 86400.) * 86400.
        else:
            dpcentdt = gr.ddt(precip_sample.p_cent) * 86400.
    
        dpcentdt_ppos = dpcentdt.where(precip_sample.p_cent>=0.)   # Find precip centroid rate where precip centroid is in the northern hemisphere
        dpcentdt_max = dpcentdt_ppos.where(dpcentdt_ppos==dpcentdt_ppos.max('xofyear'),drop=True)
        dpcentdt_max  = dpcentdt_max[0]
        pcent_max = precip_sample.p_cent.max('xofyear')
        pcent_dtmax = precip_sample.p_cent.sel(xofyear=dpcentdt_max.xofyear)    # Find the location of the preciptiation when the rate is maximum
        
        max_rate.append(float(dpcentdt_max.values))
        max_rate_lat.append(float(pcent_dtmax.values))
        max_lat.append(float(pcent_max.values))
        
    return max_rate, max_rate_lat, max_lat



if __name__ == "__main__":
    
    #data = get_data('sn_1.000')
    #max_rate, max_rate_lat, max_lat = bootstrap_method(data)
    #np.save('sn_1_bootstrap',np.array([max_rate,max_rate_lat,max_lat]))
    
    #data = get_data('sn_0.500', period_fac=0.5, months=[121,300])
    #max_rate, max_rate_lat, max_lat = bootstrap_method(data)
    #np.save('sn_p5_bootstrap',np.array([max_rate,max_rate_lat,max_lat]))
    
    #data = get_data('rt_0.500')
    #max_rate, max_rate_lat, max_lat = bootstrap_method(data)
    #np.save('rt_p5_bootstrap',np.array([max_rate,max_rate_lat,max_lat]))
    
    #data = get_data('rt_0.750')
    #max_rate, max_rate_lat, max_lat = bootstrap_method(data)
    #np.save('rt_p75_bootstrap',np.array([max_rate,max_rate_lat,max_lat]))
    
    #data = get_data('rt_1.250')
    #max_rate, max_rate_lat, max_lat = bootstrap_method(data)
    #np.save('rt_1p25_bootstrap',np.array([max_rate,max_rate_lat,max_lat]))
    
    data = get_data('rt_1.500')
    max_rate, max_rate_lat, max_lat = bootstrap_method(data)
    np.save('rt_1p5_bootstrap',np.array([max_rate,max_rate_lat,max_lat]))
    
    data = get_data('rt_1.750')
    max_rate, max_rate_lat, max_lat = bootstrap_method(data)
    np.save('rt_1p75_bootstrap',np.array([max_rate,max_rate_lat,max_lat]))
    
    data = get_data('rt_2.000')
    max_rate, max_rate_lat, max_lat = bootstrap_method(data)
    np.save('rt_2_bootstrap',np.array([max_rate,max_rate_lat,max_lat]))

