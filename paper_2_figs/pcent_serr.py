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
from pcent_rate_max import p_cent_rate_max
    
    # I want this to:
    # Load up all the years for a run
    # Rearrange them into year and pentad of year
    # Calculate the time and hemispheric mean precipitation centroid's max rate, lat of max rate, and max lat
    # Calculate the precipitation centroid's max rate, lat of max rate, and max lat for every year and both hemispheres
    # Firstly average this and compare with the value for the time and hemispheric mean
    # Then take the standard deviation of the above quantities (decide which mean to use at this point) 

def calc_standard_deviation(data, average):
    
    n = data.size
    sample_variance = ((data - average)**2.).sum(('hemisphere','year_no'))/(n-1.)
    standard_dev = np.sqrt(sample_variance.values)
    standard_error = standard_dev/np.sqrt(n)
    
    print(average, standard_dev)
    
    return standard_dev

def p_cent_error(run):
    
    max_rate_av, max_rate_lat_av, max_lat_av = p_cent_rate_max('sn_1.000')
    
    max_rate, max_rate_lat, max_lat = p_cent_rate_max(run, do_make_sym=False,months=[121,481])
    
    max_rate_mean = max_rate.mean(('hemisphere','year_no'))
    max_rate_lat_mean = max_rate_lat.mean(('hemisphere','year_no'))
    max_lat_mean = max_lat.mean(('hemisphere','year_no'))
    
    max_rate_se = calc_standard_error(max_rate, max_rate_mean.values)
    max_rate_lat_se = calc_standard_error(max_rate_lat, max_rate_lat_mean.values)
    max_lat_se = calc_standard_error(max_lat, max_lat_mean.values)
        
    max_rate_se = calc_standard_error(max_rate, max_rate_av.values)
    max_rate_lat_se = calc_standard_error(max_rate_lat, max_rate_lat_av.values)
    max_lat_se = calc_standard_error(max_lat, max_lat_av.values)
    
    #print (max_rate_av.values, max_rate.mean(('hemisphere','year_no')).values)
    #print (max_rate_lat_av.values, max_rate_lat.mean(('hemisphere','year_no')).values)
    #print (max_lat_av.values, max_lat.mean(('hemisphere','year_no')).values)
    

if __name__ == "__main__":
    
    p_cent_error('sn_1.000')