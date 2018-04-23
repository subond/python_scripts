""" 16/02/2018 Take area averages of GPCP precip over:
India: 70-100E, 5-30N
East Asia: 100-140E, 20-40N
Western North Pacific: 100-170E, 5-20N
America: 240-280E, 5-30N
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams


name_temp = '/scratch/rg419/obs_and_reanalysis/datafiles/gpcp_1dd_v1.2_p1d.%04d%02d.nc'
names = [name_temp % (m,n) for m in range( 1997, 2015) for n in range(1,13) ]

#read data into xarray 
data = xr.open_mfdataset( names, chunks={'time': 30})

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
    
    coslat = np.cos(data.precip.lat * np.pi/180.)
    sinlat = np.sin(data.precip.lat * np.pi/180.)
    
    precip_area_weighted = data.precip * coslat
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= latin[0] and data.lat[i] < latin[1]]
    lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    
    precip_mean = precip_area_weighted.sel(lon=lons).sel(lat=lats).sum(('lat','lon')) / (coslat.sel(lat=lats).sum('lat') * len(lons))

    return precip_mean


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
    



india = area_mean(region='India')
#easia = area_mean(region='EAsia')
#wnp = area_mean(region='WNP')
#america = area_mean(region='America')

india = pentad_means_of_year(india, 1997)
india.plot()

#india[0:1080].plot()
#plt.figure(2)
#easia[0:1080].plot()
#wnp[0:1080].plot()
#america[0:1080].plot()
plt.show()


