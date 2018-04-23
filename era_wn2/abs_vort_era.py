"""16/02/2018 Get absolute vorticity time series of different regions over the Asian continent... Anything interesting here?
"""

#from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
#from physics import model_constants as mc, gradients as gr
from windspharm.xarray import VectorWind
import pandas as pd


filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_vo/era_vo_200.nc'
data = xr.open_dataset(filename, chunks={'latitude': 100, 'longitude': 100})

data.load()

vo = data.vo.resample('D', dim='time', how='mean')

omega = 7.2921150e-5
f = 2 * omega * np.sin(data.latitude *np.pi/180)
    
vo = vo + f

def area_mean(region=None, lonin=[0.,360.], latin=[-90.,90.]):
        
    if region == 'India':
        lonin = [70,100]
        latin = [5,30]
    elif region == 'EAsia':
        lonin = [100,140]
        latin = [20,40]
    elif region == 'WNP':
        lonin = [100,170]
        latin = [5,20]
    elif region == 'America':
        lonin = [-120,-80]
        latin = [5,30]
    
    coslat = np.cos(vo.latitude * np.pi/180.)
    sinlat = np.sin(vo.latitude * np.pi/180.)
    
    vo_area_weighted = vo * coslat
    
    lats = [vo.latitude[i] for i in range(len(vo.latitude)) if vo.latitude[i] >= latin[0] and vo.latitude[i] < latin[1]]
    lons = [vo.longitude[i] for i in range(len(vo.longitude)) if vo.longitude[i] >= lonin[0] and vo.longitude[i] < lonin[1]]
    
    vo_mean = vo_area_weighted.sel(longitude=lons).sel(latitude=lats).sum(('latitude','longitude')) / (coslat.sel(latitude=lats).sum('latitude') * len(lons))
    
    return vo_mean


def pentad_means_of_year(data, year):
    
    data = data.sel(time=str(year))
    
    if len(data.time)==366:
        pentad = np.repeat(np.arange(1., 74.), 5)
        pentad = np.insert(pentad, 10, 2)    
    else:
        pentad = np.repeat(np.arange(1., 74.), 5)
        
    data = data.assign_coords(pentad = ('time', pentad))
    
    data = data.groupby('pentad').mean(('time'))
    
    return data
    

vo_mean = area_mean(latin=[30,45], lonin=[110,150])
vo_1997 = pentad_means_of_year(vo_mean, 1997)

plt.figure(1)
vo_1997.plot()

vo_mean = area_mean(latin=[30,45], lonin=[60,90])
vo_1997 = pentad_means_of_year(vo_mean, 1997)

plt.figure(2)
vo_1997.plot()
#plt.figure(3)
#vo_mean = area_mean(latin=[30,45], lonin=[60,90])
#vo_mean[0:1080].plot()


plt.show()