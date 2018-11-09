''' 
08/10/2018 Plot walker cell strength versus: Hadley cell strength, local Hadley cell boundary latitude, local precipitation centroid, temperature contrast
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import walker_strength, get_edge_psi

def get_lons(data, lonin):
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons

def area_mean(data, lonin, latin=[0.,90.]):
    # Take area mean over a given region
    
    lats = [data.lat[i].values for i in range(len(data.lat)) if data.lat[i] >= latin[0] and data.lat[i] <= latin[1]]
    lons = get_lons(data, lonin)
    
    coslat = np.cos(data.lat * np.pi/180.)
    
    data_weighted = data * coslat
    
    data_mean = data_weighted.sel(lon=lons).sel(lat=lats).sum(('lat','lon')) / (coslat.sel(lat=lats).sum('lat') * len(lons))
    
    return data_mean
    

def walker_strength_plot(run, cn='C0', cp='C1'):
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    walker_mag_land = walker_strength(data, lonin=[90.,180.])
    walker_mag_ocean = walker_strength(data, lonin=[270.,360.], psi_min=False)
    plt.plot(walker_mag_land.xofyear, -1.*walker_mag_land, 'x', color=cn)
    plt.plot(walker_mag_ocean.xofyear, walker_mag_ocean, 'x', color=cp)


def hadley_strength_plot(run, lonin=[370.,10.], c='C0'):
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    edge_loc, hadley_mag, psi_max_loc = get_edge_psi(data, lonin=lonin)
    plt.plot(hadley_mag.xofyear, hadley_mag, 's', color=c)

def t_surf_mean_plot(run, lonin=[-1.,361.], latin=[0.,90.], c='C0'):
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    t_surf_mean = area_mean(data.t_surf, lonin, latin=latin)
    plt.plot(t_surf_mean.xofyear, t_surf_mean, 'o', color=c)
    
def t_surf_diff_plot(run, lonin_land=[0.,180.], lonin_sea=[180.,360.], c='C0'):
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    t_surf_land = area_mean(data.t_surf, lonin_land)
    t_surf_sea = area_mean(data.t_surf, lonin_sea)
    t_surf_diff = t_surf_land - t_surf_sea
    plt.plot(t_surf_diff.xofyear, t_surf_diff, 'o', color=cd)


def walker_vs_hadley(run, lonin=[370.,10.]):
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    #walker_mag = walker_strength(data, lonin=[90.,180.])
    walker_mag = walker_strength(data, lonin=[270.,360.], psi_min=False)
    edge_loc, psi_max, psi_max_loc = get_edge_psi(data, lonin=lonin)
    times = np.intersect1d(walker_mag.xofyear.values, psi_max.xofyear)
    plt.plot(psi_max.sel(xofyear = times), walker_mag.sel(xofyear = times), 'x')



def walker_vs_temp(run, lonin=[370.,10.], lonin_land=[0.,180.], lonin_sea=[180.,360.]):
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    walker_mag = walker_strength(data, lonin=[90.,180.])
    #walker_mag = walker_strength(data, lonin=[270.,360.], psi_min=False)
    t_surf_land = area_mean(data.t_surf, lonin_land)
    t_surf_sea = area_mean(data.t_surf, lonin_sea)
    t_surf_diff = t_surf_land - t_surf_sea
    plt.figure(1)
    plt.plot(t_surf_diff.sel(xofyear=walker_mag.xofyear.values))
    plt.figure(2)
    plt.plot(walker_mag)
    plt.figure(3)
    plt.plot(t_surf_diff.sel(xofyear=walker_mag.xofyear.values), walker_mag, 'x')    


#walker_vs_hadley('half_shallow', lonin=[80.,100.])
#walker_vs_hadley('half_shallow', lonin=[170.,190.])
#walker_vs_hadley('half_shallow', lonin=[260.,280.])
#walker_vs_hadley('half_shallow', lonin=[370.,10.])

walker_strength_plot('half_shallow')
plt.figure(2)
hadley_strength_plot('half_shallow')
hadley_strength_plot('half_shallow',lonin=[170.,190.], c='C1')
hadley_strength_plot('half_shallow',lonin=[80.,100.], c='C2')
hadley_strength_plot('half_shallow',lonin=[260.,280.], c='C3')
plt.figure(3)
t_surf_mean_plot('half_shallow', lonin=[0.,180.])
t_surf_mean_plot('half_shallow', lonin=[180.,360.], c='C1')
t_surf_mean_plot('half_shallow', lonin=[0.,180.], latin=[-90.,0.])
t_surf_mean_plot('half_shallow', lonin=[180.,360.], latin=[-90.,0.], c='C1')
t_surf_mean_plot('ap_2', c='C2', latin=[-90.,90.])
t_surf_mean_plot('ap_20', c='C3', latin=[-90.,90.])
plt.show()