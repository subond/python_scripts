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


def get_data(run, months=[121,481], period_fac=1., days=False, filename = 'plev_pentad'):
    # Load in the data, reshape it, flip it and add flipped data onto the end
    # This way, all means will be hemispheric and time, and both hemispheres will be included in sampling
    print(run)
    data = precip_load_and_reshape(run, months=months, period_fac=period_fac, days=days, filename=filename)
    
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
    eq_rate = []
    
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
        
        p_cent_abs = np.abs(precip_sample.p_cent.where(dpcentdt>=0.))        
        dpdt_eq = (dpcentdt.where(p_cent_abs == p_cent_abs.min('xofyear'), drop=True)).values[0]
        #dpdt_eq = dpcentdt.where(precip_sample.p_cent == precip_sample.p_cent.min('xofyear'), drop=True).max()
        
        max_rate.append(float(dpcentdt_max.values))
        max_rate_lat.append(float(pcent_dtmax.values))
        max_lat.append(float(pcent_max.values))
        eq_rate.append(float(dpdt_eq))
        
    return max_rate, max_rate_lat, max_lat, eq_rate



if __name__ == "__main__":
    
    
    #data = get_data('sn_0.500', period_fac=0.5, months=[121,300])
    #max_rate, max_rate_lat, max_lat, eq_rate = bootstrap_method(data)
    #np.save('sn_p5_bootstrap',np.array([max_rate,max_rate_lat,max_lat, eq_rate]))
    
    #runs = ['mld_2.5', 'mld_5', 'sn_1.000', 'mld_15', 'mld_20']
            
    #for run in runs:
    #    try:
    #        data = get_data(run)
    #        max_rate, max_rate_lat, max_lat, eq_rate = bootstrap_method(data)
    #        np.save(run + '_bootstrap',np.array([max_rate,max_rate_lat,max_lat, eq_rate]))
    #        data.close()
    #    except:
    #        print(run + ' failed')
    
    
    #data = get_data('sn_0.500', period_fac=0.5, months=[121,301])
    #max_rate, max_rate_lat, max_lat, eq_rate = bootstrap_method(data)
    #np.save('sn_0.500_bootstrap',np.array([max_rate,max_rate_lat,max_lat, eq_rate]))
    
    #data = get_data('sn_2.000', period_fac=2., months=[121,841])
    #max_rate, max_rate_lat, max_lat, eq_rate = bootstrap_method(data)
    #np.save('sn_2.000_bootstrap',np.array([max_rate,max_rate_lat,max_lat, eq_rate]))
    
    #data = get_data('sn_0.250', period_fac=0.25, months=[121,211], days=True, filename='plev_daily_mean')
    #max_rate, max_rate_lat, max_lat, eq_rate = bootstrap_method(data, days=True)
    #np.save('sn_0.250_bootstrap',np.array([max_rate,max_rate_lat,max_lat, eq_rate]))
    
    #data = get_data('sn_3.000', period_fac=3., months=[145,1225])
    #max_rate, max_rate_lat, max_lat, eq_rate = bootstrap_method(data)
    #np.save('sn_3.000_bootstrap',np.array([max_rate,max_rate_lat,max_lat, eq_rate]))
    
    data = get_data('sn_4.000', period_fac=4., months=[145,1585])
    max_rate, max_rate_lat, max_lat, eq_rate = bootstrap_method(data)
    np.save('sn_4.000_bootstrap',np.array([max_rate,max_rate_lat,max_lat, eq_rate]))
    
