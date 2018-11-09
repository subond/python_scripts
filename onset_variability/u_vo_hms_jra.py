'''04/10/2018 produce yearly u and vo 200 hPa hovmoller plots for JRA-55 data, with onset date marked on'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm
from data_handling_updates import gradients as gr
from windspharm.xarray import VectorWind

plot_dir = '/scratch/rg419/plots/onset_variability/u_vo_hms_jra/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

# Load in onset dates
onsets_scsm = np.load('/scratch/rg419/python_scripts/onset_variability/jra_onsets_scsm.npy')
onsets_scsm = xr.DataArray(onsets_scsm, coords={'year': ('year', range(1958,2017))}, dims=['year'])

years = range(1958,2017)

# Load in JRA data
#data_u = xr.open_dataset('/disca/share/rg419/jra_ucomp_daily_200.nc')
data_u = xr.open_dataset('/disca/share/rg419/jra_ucomp_daily_850.nc')
data_u = data_u['var33'].load().loc['1958-01':'2016-12']
#data_u_clim = data_u.groupby('time.dayofyear').mean('time')
#print(data_u_clim)

# v has different time coord to u, presumably due to how Stephen has downloaded/averaged. I think the two are equivalent, so just substitute the time dimension into v
data_v_temp = xr.open_dataset('/disca/share/rg419/jra_vcomp_daily_200.nc')
data_v = xr.DataArray(data_v_temp.var34.values, coords=data_u.coords, dims=data_u.dims)

dvdx = gr.ddx(data_v, a=6371.0e3) 
dudy = gr.ddy(data_u, a=6371.0e3)
data_vo = (dvdx - dudy)*86400.

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

def make_pentad_clim(data,years):
    data_temp = np.zeros([len(years),73,73,144])
    for i in range(len(years)):
        data_temp[i,:,:,:] = pentad_means_of_year(data, years[i]).values
    
    data_temp = xr.DataArray(data_temp, coords = {'year': years, 'pentad': np.arange(1., 74.), 'lat': data_u.lat, 'lon': data_u.lon}, 
                          dims=['year','pentad','lat','lon'])
    data_clim = data_temp.mean('year')
    return data_clim
   
        
def plot_hm(data, title, levels=np.arange(-15.,16.,1.), lonin=[100,140], anom=True):
        
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    # Calculate pentad climatology
    data_clim = make_pentad_clim(data,years)
    
    data = data.sel(lon=lons).mean('lon')
    
    # Start figure with 60 subplots
    rcParams['figure.figsize'] = 30, 20
    rcParams['font.size'] = 14
    fig, axes = plt.subplots(6, 10, sharex='col', sharey='row')
    axlist = [element for tupl in axes for element in tupl]
    
    for i in range(len(years)):
        if anom:
            data_pentad = pentad_means_of_year(data, years[i]) - data_clim.sel(lon=lons).mean('lon')
        else:
            data_pentad = pentad_means_of_year(data, years[i])
        f1 = data_pentad.plot.contourf(ax=axlist[i], x='pentad',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
        axlist[i].plot([onsets_scsm[i], onsets_scsm[i]], [-60,60], 'k')
        axlist[i].set_ylim(0.,60.)
        axlist[i].set_xlim(0.,40.)
        axlist[i].set_yticks(np.arange(0.,61.,10.))
        axlist[i].set_title(str(years[i]))
        axlist[i].grid(True,linestyle=':')

        
    plt.subplots_adjust(left=0.02, right=0.99, top=0.98, bottom=0.05, hspace=0.3, wspace=0.2)
    
    cb1=fig.colorbar(f1, ax=axlist, use_gridspec=True, orientation = 'horizontal',fraction=0.02, pad=0.05, aspect=60, shrink=0.75)
    
    plt.savefig(plot_dir + title + '.pdf', format='pdf')
    plt.close()


#plot_hm(data_u, 'ucomp_hm')
#plot_hm(data_vo, 'vo_hm', levels=np.arange(-2.,2.1,0.2))
#plot_hm(data_u, 'ucomp_full_hm', anom=False, levels=np.arange(-40.,41.,5.))
plot_hm(data_u, 'ucomp_full_hm_850', anom=False, levels=np.arange(-14.,15.,2.))
#plot_hm(data_vo, 'vo_full_hm', anom=False, levels=np.arange(-5.,5.1,0.5))

data_u.close()
data_v.close()
data_vo.close()



