'''28/09/2018 Calculate onset dates from daily NCEP-NCAR reanalysis, and save'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
import pandas as pd
from climatology import bob_onset, scsm_onset, ism_onset


filename = '/disca/share/rg419/ncep_ncar_uwnd_daily_850.nc'   # Have taken 850hPa level from Stephen's datafile and saved seperately 

data = xr.open_dataset('/disca/share/rg419/ncep_ncar_uwnd_daily_850_1948_2016.nc')

#data_1 = xr.open_dataset('/disca/share/rg419/ncep_ncar_uwnd_daily_850_1948_1987.nc')
#data_2 = xr.open_dataset('/disca/share/rg419/ncep_ncar_uwnd_daily_850_1988_2017.nc')
#data = xr.open_mfdataset(['/disca/share/rg419/ncep_ncar_uwnd_daily_850_1948_1987.nc',
#                          '/disca/share/rg419/ncep_ncar_uwnd_daily_850_1988_2017.nc'])
#data = xr.concat([data_1, data_2], dim='time')

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
onset_ncep=[]    
for year in range(1948,2017):
    print(year)
    u_year = pentad_means_of_year(data.uwnd, year)
    u_ncep = scsm_onset(u_year, pdim='lev', pentaddim='pentad', latdim='lat', londim='lon')
    onset_ncep.append(u_ncep[1])

print(onset_ncep)
np.save('ncep_onsets_scsm',np.array(onset_ncep))


# Get BoB onset date
onset_ncep=[]    
for year in range(1948,2017):
    print(year)
    u_year = day_means_of_year(data.uwnd, year)
    u_ncep = bob_onset(u_year, pdim='lev', daydim='day', latdim='lat', londim='lon')
    onset_ncep.append(u_ncep[1])

print(onset_ncep)
np.save('ncep_onsets_bob_days',np.array(onset_ncep))


# Get ISM onset date
onset_ncep=[]    
for year in range(1948,2017):
    print(year)
    u_year = day_means_of_year(data.uwnd, year)
    u_ncep = ism_onset(u_year, pdim='lev', daydim='day', latdim='lat', londim='lon')
    onset_ncep.append(u_ncep[1])

print(onset_ncep)
np.save('ncep_onsets_ism_days',np.array(onset_ncep))

