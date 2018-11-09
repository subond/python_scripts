'''01/10/2018 Calculate (updated) onset dates from daily ERA-Interim reanalysis, and save'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
import pandas as pd
from climatology import bob_onset, scsm_onset, ism_onset


data = xr.open_dataset('/disca/share/rg419/era_ucomp_daily_850.nc')


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

# Get SCSM onset date 
onset_era=[]    
for year in range(1979,2017):
    print(year)
    u_year = pentad_means_of_year(data.u, year)
    u_era = scsm_onset(u_year, pdim='level', pentaddim='pentad', latdim='latitude', londim='longitude')
    onset_era.append(u_era[1])

print(onset_era)
np.save('era_onsets_scsm',np.array(onset_era))


# Get BoB onset date
onset_era=[]    
for year in range(1979,2017):
    print(year)
    u_year = day_means_of_year(data.u, year)
    u_era = bob_onset(u_year, pdim='level', daydim='day', latdim='latitude', londim='longitude')
    onset_era.append(u_era[1])

print(onset_era)
np.save('era_onsets_bob_days',np.array(onset_era))


# Get ISM onset date
onset_era=[]    
for year in range(1979,2017):
    print(year)
    u_year = day_means_of_year(data.u, year)
    u_era = ism_onset(u_year, pdim='level', daydim='day', latdim='latitude', londim='longitude')
    onset_era.append(u_era[1])

print(onset_era)
np.save('era_onsets_ism_days',np.array(onset_era))

