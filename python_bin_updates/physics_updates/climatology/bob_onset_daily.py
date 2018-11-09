""" 2/08/2018 Function to calculate the Mao and Wu 2006 Bay of Bengal monsoon onset diagnostic.
Takes in u field, and averages this between 5 and 15N, and 90-100 E if a lon dimension is given. 
Time must be in days
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from data_handling_updates import isca_load_and_reshape


def bob_onset(u, daydim='xofyear', pdim='pfull', latdim='lat', londim='lon', print_onset=False, lev=850.):
    
    # Check if the input day name is in the xarray dimensions, rename to day if so
    if daydim in u.dims:
        u = u.rename(new_name_or_name_dict={daydim: 'day'})
    else:    # Otherwise cause a crash 
        u['daydim']
    if latdim in u.dims: # Same with lat
        u = u.rename(new_name_or_name_dict={latdim: 'lat'})
    else:    
        u['latdim']
    if londim in u.dims: # Same with lon
        u = u.rename(new_name_or_name_dict={londim: 'lon'})
    else:    
        u['londim']
        
    # Select 850 hPa pressure level, if multiple pressure levels are present
    if pdim in u.dims:
        u = u.rename(new_name_or_name_dict={pdim: 'pfull'})
        u = u.sel( pfull = lev)
    
    coslat = np.cos(u.lat * np.pi/180.)
    sinlat = np.sin(u.lat * np.pi/180.)
    
    u_area_weighted = u * coslat # Weight u by coslat, to do area weighting
    
    # Get specified lats and lons, select day range to look at
    lats = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= 5. and u.lat[i] <= 15.]
    lons = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= 90. and u.lon[i] <= 100.]
    days = [u.day[i].values for i in range(len(u.day)) if u.day[i] >=  90]
        
    # Calculate the area mean u
    u_mean = u_area_weighted.sel(lon=lons).sel(lat=lats).sum(('lat','lon')) / (coslat.sel(lat=lats).sum('lat') * len(lons))
        
    for i in range(len(days)):  # For each day
        if (all(u_mean.sel(day=days[i:i+11]).values > 0)
            and (u_mean.sel(day=days[i-1]).values < 0)): # Check if u_mean changes from negative to positive, and if it stays so for the next 10 days
            onset_day = days[i]  # If all that is true, that's your onset day
            if print_onset:
                print('Onset day: ', onset_day) # Print it
            return u_mean, onset_day # Return it, and u
    if print_onset:
        print('No onset') # If you loop the whole way through and there's no onset, i.e. code still hasn't returned an onset date and finished, print that if wanted
    return u_mean, None # And return u and an empty value

    

if __name__ == "__main__":
    
    filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_850.nc'
    data = xr.open_dataset(filename, chunks={'latitude': 100, 'longitude': 100})
    data = data.resample('D', dim='time', how='mean')
        
    def day_means_of_year(data, year):
        data_year = data.sel(time=str(year))
        if len(data_year.time)==366:
            day = np.arange(1.,367.)   
        else:
            day = np.arange(1.,366.)
        data_year = data_year.assign_coords(day = ('time', day))
        data_year = data_year.groupby('day').mean(('time'))
        return data_year
        
    onset_era=[]    
    for year in range(1979,2017):
        print(year)
        u_year = day_means_of_year(data.u, year)
        u_era = bob_onset(u_year, pdim='level', daydim='day', latdim='latitude', londim='longitude')
        onset_era.append(u_era[1])
    
    print(onset_era)
    np.save('era_onsets_bob_days',np.array(onset_era))
    