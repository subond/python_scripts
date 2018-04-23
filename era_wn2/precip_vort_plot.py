""" 16/02/2018 Plot area averaged precip over India and East Asia on same time axis as area average 200 hPa absolute vorticity to the north of this region
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
import pandas as pd

# Load datasets
name_temp = '/scratch/rg419/obs_and_reanalysis/datafiles/gpcp_1dd_v1.2_p1d.%04d%02d.nc'
names = [name_temp % (m,n) for m in range( 1997, 2015) for n in range(1,13) ]
data_precip = xr.open_mfdataset( names, chunks={'time': 30})

filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_vo/era_vo_200.nc'
data_vo = xr.open_dataset(filename, chunks={'latitude': 100, 'longitude': 100})
data_vo.load()

# Get daily mean absolute vorticity
vo = data_vo.vo.resample('D', dim='time', how='mean')
omega = 7.2921150e-5
f = 2 * omega * np.sin(data_vo.latitude *np.pi/180)
vo = vo + f
vo = vo.sel(level=200.) * 86400.


# Define useful functions - area mean and pentad means

def area_mean(data, region=None, lonin=[0.,360.], latin=[-90.,90.]):
    
    if 'latitude' in data.dims:
        data = data.rename({'latitude': 'lat'})
    if 'longitude' in data.dims:
        data = data.rename({'longitude': 'lon'})
    
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
    
    coslat = np.cos(data.lat * np.pi/180.)
    sinlat = np.sin(data.lat * np.pi/180.)
    
    data_area_weighted = data * coslat
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= latin[0] and data.lat[i] < latin[1]]
    lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    
    data_mean = data_area_weighted.sel(lon=lons).sel(lat=lats).sum(('lat','lon')) / (coslat.sel(lat=lats).sum('lat') * len(lons))

    return data_mean


def pentad_means_of_year(data, year):
    
    data = data.sel(time=str(year))
    
    if len(data.time)==366:
        pentad = np.repeat(np.arange(1., 74.), 5)
        pentad = np.insert(pentad, 10, 2)    
    else:
        pentad = np.repeat(np.arange(1., 74.), 5)
        
    data = data.assign_coords(pentad = ('time', pentad))
    
    data = data.groupby('pentad').mean(('time'))
        
    return data.values
    

def pentad_datetime(data, years):
    # Set up pentads date time object
    pentads =[]
    for year in years:
        for i in range(73):
            date = str(data.time[i*5].values)[4:] 
            pentads.append(str(year) + date) 
    pentads = pd.to_datetime(pentads)
    return pentads
    

def get_pentads(data, years):
    
    data_out = np.zeros(len(years)*73,)
    i=0
    for year in years:
        data_out[i*73: (i+1)*73] = pentad_means_of_year(data, year)
        i=i+1
    
    pentads = pentad_datetime(data, years)
    data_out = xr.DataArray(data_out, coords=[pentads], dims=['time'])
    
    return data_out
    
    
    
precip_mean_india = area_mean(data_precip.precip, region='India')
precip_mean_easia = area_mean(data_precip.precip, region='EAsia')

precip_india = get_pentads(precip_mean_india, range(1997,2015))
precip_easia = get_pentads(precip_mean_easia, range(1997,2015))

vo_mean = area_mean(vo, latin=[30,45], lonin=[110,150])
vo_E = get_pentads(vo_mean, range(1997,2015))

vo_mean = area_mean(vo, latin=[30,45], lonin=[60,90])
vo_W = get_pentads(vo_mean, range(1997,2015))

precip_clim_india = precip_india.groupby('time.dayofyear').mean('time')
precip_clim_easia = precip_easia.groupby('time.dayofyear').mean('time')
vo_clim_E = vo_E.groupby('time.dayofyear').mean('time')
vo_clim_W = vo_W.groupby('time.dayofyear').mean('time')


#print vo_E.time['month']

rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 14
plot_dir = '/scratch/rg419/plots/era_wn2/'

#f, ax = plt.subplots(1, 1)
fig, (ax1, ax2) = plt.subplots(2)
precip_india[0:657].plot(ax=ax1)
ax_twin = ax1.twinx()
vo_E[0:657].plot(ax=ax_twin, color='k')
precip_india[657:].plot(ax=ax2)
ax2_twin = ax2.twinx()
vo_E[657:].plot(ax=ax2_twin, color='k')
ax1.set_ylabel('Precip, mm/day')
ax2.set_ylabel('Precip, mm/day')
ax_twin.set_ylabel('Absolute vorticity, /day')
ax2_twin.set_ylabel('Absolute vorticity, /day')
ax1.spines['left'].set_color('blue')
ax2.spines['left'].set_color('blue')
ax1.set_title('Indian precip, eastern vorticity')
plt.savefig(plot_dir + 'indian_precip_E_vor.pdf', format='pdf')
plt.close()

fig, (ax1, ax2) = plt.subplots(2)
precip_india[0:657].plot(ax=ax1)
ax_twin = ax1.twinx()
vo_W[0:657].plot(ax=ax_twin, color='k')
precip_india[657:].plot(ax=ax2)
ax2_twin = ax2.twinx()
vo_W[657:].plot(ax=ax2_twin, color='k')
ax1.set_ylabel('Precip, mm/day')
ax2.set_ylabel('Precip, mm/day')
ax_twin.set_ylabel('Absolute vorticity, /day')
ax2_twin.set_ylabel('Absolute vorticity, /day')
ax1.spines['left'].set_color('blue')
ax2.spines['left'].set_color('blue')
ax1.set_title('Indian precip, western vorticity')
plt.savefig(plot_dir + 'indian_precip_W_vor.pdf', format='pdf')
plt.close()

fig, (ax1, ax2) = plt.subplots(2)
precip_easia[0:657].plot(ax=ax1)
ax_twin = ax1.twinx()
vo_E[0:657].plot(ax=ax_twin, color='k')
precip_easia[657:].plot(ax=ax2)
ax2_twin = ax2.twinx()
vo_E[657:].plot(ax=ax2_twin, color='k')
ax1.set_ylabel('Precip, mm/day')
ax2.set_ylabel('Precip, mm/day')
ax_twin.set_ylabel('Absolute vorticity, /day')
ax2_twin.set_ylabel('Absolute vorticity, /day')
ax1.spines['left'].set_color('blue')
ax2.spines['left'].set_color('blue')
ax1.set_title('East Asian precip, eastern vorticity')
plt.savefig(plot_dir + 'easian_precip_E_vor.pdf', format='pdf')
plt.close()

fig, (ax1, ax2) = plt.subplots(2)
precip_easia[0:657].plot(ax=ax1)
ax_twin = ax1.twinx()
vo_W[0:657].plot(ax=ax_twin, color='k')
precip_easia[657:].plot(ax=ax2)
ax2_twin = ax2.twinx()
vo_W[657:].plot(ax=ax2_twin, color='k')
ax1.set_ylabel('Precip, mm/day')
ax2.set_ylabel('Precip, mm/day')
ax_twin.set_ylabel('Absolute vorticity, /day')
ax2_twin.set_ylabel('Absolute vorticity, /day')
ax1.spines['left'].set_color('blue')
ax2.spines['left'].set_color('blue')
ax1.set_title('East Asian precip, western vorticity')
plt.savefig(plot_dir + 'easian_precip_W_vor.pdf', format='pdf')
plt.close()


rcParams['figure.figsize'] = 8, 5
fig, ax = plt.subplots(1,1)
precip_clim_easia.plot(ax=ax)
ax_twin = ax.twinx()
vo_clim_E.plot(ax=ax_twin, color='k')
ax.set_ylabel('Precip, mm/day')
ax_twin.set_ylabel('Absolute vorticity, /day')
ax.set_ylim([0,10])
ax_twin.set_ylim([5,11])
ax.spines['left'].set_color('blue')
ax.set_title('East Asian precip, eastern vorticity')
plt.savefig(plot_dir + 'easian_precip_E_vor_clim.pdf', format='pdf')
plt.close()

fig, ax = plt.subplots(1,1)
precip_clim_easia.plot(ax=ax)
ax_twin = ax.twinx()
vo_clim_W.plot(ax=ax_twin, color='k')
ax.set_ylabel('Precip, mm/day')
ax_twin.set_ylabel('Absolute vorticity, /day')
ax.set_ylim([0,10])
ax_twin.set_ylim([5,11])
ax.spines['left'].set_color('blue')
ax.set_title('East Asian precip, western vorticity')
plt.savefig(plot_dir + 'easian_precip_W_vor_clim.pdf', format='pdf')
plt.close()

fig, ax = plt.subplots(1,1)
precip_clim_india.plot(ax=ax)
ax_twin = ax.twinx()
vo_clim_E.plot(ax=ax_twin, color='k')
ax.set_ylabel('Precip, mm/day')
ax_twin.set_ylabel('Absolute vorticity, /day')
ax.set_ylim([0,10])
ax_twin.set_ylim([5,11])
ax.spines['left'].set_color('blue')
ax.set_title('Indian precip, eastern vorticity')
plt.savefig(plot_dir + 'indian_precip_E_vor_clim.pdf', format='pdf')
plt.close()

fig, ax = plt.subplots(1,1)
precip_clim_india.plot(ax=ax)
ax_twin = ax.twinx()
vo_clim_W.plot(ax=ax_twin, color='k')
ax.set_ylabel('Precip, mm/day')
ax_twin.set_ylabel('Absolute vorticity, /day')
ax.set_ylim([0,10])
ax_twin.set_ylim([5,11])
ax.spines['left'].set_color('blue')
ax.set_title('Indian precip, western vorticity')
plt.savefig(plot_dir + 'indian_precip_W_vor_clim.pdf', format='pdf')
plt.close()



