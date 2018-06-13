"""
Calculate rate of movement of precipitation centroid for different seasons 02/02/2018

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_centroid
from data_handling_updates import gradients as gr
from pylab import rcParams
from data_handling_updates import make_sym


def precip_load_and_reshape(run, months=[121,181], period_fac=1., filename='plev_pentad', days=False):
    '''Load and reshape multi-year data into year/day of year shape'''
    name_temp = '/disca/share/rg419/Data_moist/' + run + '/run%04d/' + filename + '.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    
    #read data into xarray 
    data = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    
    try:
        data = data.condensation_rain + data.convection_rain
    except:
        data = data.precipitation
    
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
    
    data_rs = xr.Dataset({'precipitation': data_rs}, coords=data_rs.coords)
    data_rs = data_rs.transpose('year_no', 'xofyear', 'lat', 'lon')

    return data_rs
    

def p_cent_rate_max(runs, do_make_sym=True, months=None, days=None, period_fac=1.):
    # Get the maximum rate of change of the precipitation centroid latitude, and the latitude at which this occurs.
    
    # Ensure runs is a list not a string
    if not isinstance(runs,list):
        runs=[runs]
    if not isinstance(period_fac,list):
        period_fac=[period_fac]*len(runs)
    # Set up empty lists
    max_rate = []
    max_rate_lat = []
    max_lat = []
    if not do_make_sym: # If we're not averaging both hemispheres, get southern peak magnitude and location too
        max_rate_s = []
        max_rate_lat_s = []
        max_lat_s = []
    
    if days==None:
        days=[False]*len(runs)
        
    j=0
    for run in runs: 
        if months==None:   # If no months provided, look for a climatology
            # Open dataset
            data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
            # Get total precip
            try:
                data['precipitation'] = data.condensation_rain + data.convection_rain
            except:
                data['precipitation'] = data.precipitation
        else: # If months are provided, load the data for these, and reshape to years/day of year
            try: # First make sure months is a list of lists
                len(months[0])
            except:
                months = [months]
            data = precip_load_and_reshape(run, months=months[j], period_fac=period_fac[j])
        if do_make_sym: # If symmetric data is wanted, average over both hemispheres (NB currently only set up for a climatology 2/05/18)
            data['precipitation'] = make_sym(data.precipitation)

        # Locate precipitation centroid
        precip_centroid(data)
        
        # Get rate of movement of precip centroid
        if days[j]:
            dpcentdt = gr.ddt(data.p_cent, secperunit = 86400.) * 86400.
        else:
            dpcentdt = gr.ddt(data.p_cent) * 86400.
                
        
        def get_max(dpcentdt_pmask, sh=False):  # Function to find the maximum ratea and it's latitude and the maximum latitude.
            dpcentdt_max = dpcentdt_pmask.where(dpcentdt_pmask==dpcentdt_pmask.max('xofyear'),drop=True)   # Find the maximum rate
            numdims = len(dpcentdt_max.shape) # Check if this is just 1 year or if it's many years
            if numdims != 1: # If it's many years, evaluate the metrics for each year
                dpdt_max_temp=[] #First set up empty lists
                pcent_dtmax=[]
                pcent_max=[]
                for y in range(dpcentdt_max.shape[0]): # Loop over years
                    if sh: #If it's the southern hemisphere, look for the min latitude of pcent, and make it positive
                        #data.p_cent.sel(year_no=y).plot()
                        #plt.show()
                        pcent_max.append(-1.*data.p_cent.sel(year_no=y).min('xofyear'))
                    else: # Otherwise look for the maximum latitude
                        pcent_max.append(data.p_cent.sel(year_no=y).max('xofyear'))
                    # Find the maximum rate for the given year, remove any nans that have been added. If there are multiple values, use the smallest
                    dpdt_max_temp.append(dpcentdt_max[y,:].dropna('xofyear')[0])
                    # Find the location of the max rate for that year
                    pcent_dtmax.append(data.p_cent.sel(year_no=y,xofyear=dpdt_max_temp[y].xofyear.values))
                # Convert all results to datarrays
                pcent_dtmax = xr.DataArray(np.asarray(pcent_dtmax), coords=[data.year_no.values], dims=['year_no'])
                dpcentdt_max = xr.DataArray(np.asarray(dpdt_max_temp), coords=[data.year_no.values], dims=['year_no'])
                pcent_max = xr.DataArray(np.asarray(pcent_max), coords=[data.year_no.values], dims=['year_no'])
            
            else: # If it's only 1D
                if len(dpcentdt_max) > 1: # If there are multiple values for the maximum take the smallest
                    dpcentdt_max  = dpcentdt_max[0].expand_dims('xofyear',axis=0) 
                if sh: # If it's the southern hemisphere find the min lat of pcent, and make it positive
                    pcent_max = -1.*data.p_cent.min('xofyear')
                else: # Otherwise find the maximum latitude
                    pcent_max = data.p_cent.max('xofyear')
                pcent_max = np.expand_dims(pcent_max, axis=1)
                pcent_dtmax = data.p_cent.sel(xofyear=dpcentdt_max.xofyear)    # Find the location of the preciptiation when the rate is maximum
            return dpcentdt_max, pcent_dtmax, pcent_max
        
        
        dpcentdt_ppos = dpcentdt.where(data.p_cent>=0.)   # Find precip centroid rate where precip centroid is in the northern hemisphere
        dpcentdt_max, pcent_dtmax, pcent_max = get_max(dpcentdt_ppos) # Get the maximum etc
        max_rate.append(dpcentdt_max) # Add values to lists for different runs
        max_rate_lat.append(pcent_dtmax)
        max_lat.append(pcent_max)
        
        if not do_make_sym:
            dpcentdt_pneg = -1.*dpcentdt.where(data.p_cent<=0.)   # Find precip centroid rate where precip centroid is in the southern hemisphere, and make poleward positive
            dpcentdt_max_s, pcent_dtmax_s, pcent_max_s = get_max(dpcentdt_pneg, sh=True) # Get the maximum etc
            max_rate_s.append(dpcentdt_max_s) # Add values to lists for different runs
            max_rate_lat_s.append(-1.*pcent_dtmax_s)
            max_lat_s.append(pcent_max_s)
        
        j=j+1 # Add 1 to counter to get values for next run
    
    # Convert all output to xarrays   
    if do_make_sym: 
        if months==None:
            coords_yrno = ['clim']
        else:
            coords_yrno = data.year_no.values
        max_rate = xr.DataArray(np.asarray(max_rate), coords=[runs,coords_yrno], dims=['run','year_no'])
        max_rate_lat = xr.DataArray(np.asarray(max_rate_lat), coords=[runs,coords_yrno], dims=['run','year_no'])
        max_lat = xr.DataArray(np.asarray(max_lat), coords=[runs,coords_yrno], dims=['run','year_no'])

    else:
        if months==None:
            coords_yrno = ['clim']
        else:
            coords_yrno = data.year_no.values
        max_rate = (xr.DataArray(np.asarray([max_rate,max_rate_s]), coords=[['n','s'],runs,coords_yrno], dims=['hemisphere','run','year_no'])).transpose('run','year_no','hemisphere')
        max_rate_lat = (xr.DataArray(np.asarray([max_rate_lat,max_rate_lat_s]), coords=[['n','s'],runs,coords_yrno], dims=['hemisphere','run','year_no'])).transpose('run','year_no','hemisphere')
        max_lat = (xr.DataArray(np.asarray([max_lat,max_lat_s]), coords=[['n','s'],runs,coords_yrno], dims=['hemisphere','run','year_no'])).transpose('run','year_no','hemisphere')

    
    return max_rate, max_rate_lat, max_lat

        

if __name__ == "__main__":
    
    max_rate_av, max_rate_lat_av, max_lat_av = p_cent_rate_max(['rt_0.500', 'rt_0.750', 'sn_1.000','rt_1.250', 'rt_1.500', 'rt_1.750'])
    max_rate, max_rate_lat, max_lat = p_cent_rate_max(['sn_1.000'], do_make_sym=False,period_fac=1., months=[121,481])
    max_rate, max_rate_lat, max_lat = p_cent_rate_max(['sn_1.000','sn_0.500'], do_make_sym=False,period_fac=[1.,0.5], months=[[121,481],[121,301]])
    
    print (max_rate_av.values, max_rate.mean(('hemisphere','year_no')).values)
    print (max_rate_lat_av.values, max_rate_lat.mean(('hemisphere','year_no')).values)
    print (max_lat_av.values, max_lat.mean(('hemisphere','year_no')).values)
    