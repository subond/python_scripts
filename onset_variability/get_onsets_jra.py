'''4/09/2018 Calculate onset dates from daily JRA reanalysis, and save'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
import pandas as pd
from climatology import bob_onset, scsm_onset, ism_onset


filename = '/disca/share/rg419/jra_ucomp_daily_850.nc'   # Have taken 850hPa level from Stephen's datafile and saved seperately 

data = xr.open_dataset(filename)
#data = data.sel(lev=85000.)
#data.to_netcdf('/disca/share/rg419/jra_ucomp_daily_850.nc')

def day_means_of_year(data, year):   # Function to get day of year
    data_year = data.sel(time=str(year))
    if len(data_year.time)==366:
        day = np.arange(1.,367.)   
    else:
        day = np.arange(1.,366.)
    data_year = data_year.assign_coords(day = ('time', day))
    data_year = data_year.groupby('day').mean(('time'))
    return data_year


def pentad_means_of_year(data, year):  # Function to get pentad of year
        data_year = data.sel(time=str(year))
        if len(data_year.time)==366:
            pentad = np.repeat(np.arange(1., 74.), 5)
            pentad = np.insert(pentad, 10, 2)    
        else:
            pentad = np.repeat(np.arange(1., 74.), 5)
        data_year = data_year.assign_coords(pentad = ('time', pentad))
        data_year = data_year.groupby('pentad').mean(('time'))
        return data_year


def pentad_mean_climatology(data, years):  # Function to get pentad of year
    pentad_years = np.array([])
    for year in years:
        data_year = data.sel(time=str(year))
        if len(data_year.time)==366:
            pentad = np.repeat(np.arange(1., 74.), 5)
            pentad = np.insert(pentad, 10, 2)    
        else:
            pentad = np.repeat(np.arange(1., 74.), 5)    
        np.concatenate((pentad_years, pentad))
        print(pentad_years.shape)
        
    data = data.assign_coords(pentad = ('time', pentad_years))
    print(data.pentad_years.values[0:720])
    
    data = data.groupby('pentad').mean(('time'))
    
    return data_pentads


# Get BoB onset date
#onset_jra=[]    
#for year in range(1958,2017):
#    print(year)
#    u_year = day_means_of_year(data.var33, year)
#    u_jra = bob_onset(u_year, pdim='lev', daydim='day', latdim='lat', londim='lon')
#    onset_jra.append(u_jra[1])

#print(onset_jra)
#np.save('jra_onsets_bob_days',np.array(onset_jra))


# Get ISM onset date
#onset_jra=[]    
#for year in range(1958,2017):
#    print(year)
#    u_year = day_means_of_year(data.var33, year)
#    u_jra = ism_onset(u_year, pdim='lev', daydim='day', latdim='lat', londim='lon')
#    onset_jra.append(u_jra[1])

#print(onset_jra)
#np.save('jra_onsets_ism_days',np.array(onset_jra))


# Get SCSM onset date 
onset_jra=[]    
for year in range(1958,2017):
    print(year)
    u_year = pentad_means_of_year(data.var33, year)
    u_jra = scsm_onset(u_year, pdim='lev', pentaddim='pentad', latdim='lat', londim='lon')
    onset_jra.append(u_jra[1])

print(onset_jra)
np.save('jra_onsets_scsm',np.array(onset_jra))