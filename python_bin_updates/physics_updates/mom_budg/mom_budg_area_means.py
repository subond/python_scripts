""" 03/04/2018 Take area averages of Isca precip over:
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
from mom_budg import mom_budg_fun


def area_mean(run, var, region=None, lonin=[0.,360.], latin=[-90.,90.]):
        
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
    
    data = mom_budg_fun(run)
    
    coslat = np.cos(data[var].lat * np.pi/180.)
    sinlat = np.sin(data[var].lat * np.pi/180.)
    
    var_area_weighted = data[var] * coslat
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= latin[0] and data.lat[i] < latin[1]]
    lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    
    var_mean = var_area_weighted.sel(lon=lons).sel(lat=lats).sum(('lat','lon')) / (coslat.sel(lat=lats).sum('lat') * len(lons))

    return var_mean


        

if __name__ == "__main__":
    
    rcParams['figure.figsize'] = 5, 3
    rcParams['font.size'] = 14
    
    plot_dir = '/scratch/rg419/plots/climatology/mom_budg_area_means/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    india_ageo_ctrl = area_mean('control_qflux', 'fv_ageo', region='India')
    india_ageo_no_TIP = area_mean('no_TIP', 'fv_ageo', region='India')
    india_geo_ctrl = area_mean('control_qflux', 'dphidx', region='India')
    india_geo_no_TIP = area_mean('no_TIP', 'dphidx', region='India')
    
    fig = plt.figure()
    ax = plt.subplot(111)
    india_ageo_ctrl.plot(color='k', linewidth=2)
    (-1.*india_geo_ctrl).plot(color='k', linestyle='--', linewidth=2)
    india_ageo_no_TIP.plot(color='r', linewidth=2)
    (-1.*india_geo_no_TIP).plot(color='r', linestyle='--', linewidth=2)
    plt.xlabel('Pentad')
    plt.ylabel('Momentum tendency, m/s/day')
    plt.xlim([1,72])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    #ax.legend(['Control', 'No Tibet'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
    plt.subplots_adjust(left=0.15, right=0.75, top=0.95, bottom=0.2)
    plt.savefig(plot_dir + 'fv_india_no_TIP.pdf', format='pdf')
    plt.close()
    
    
    easia_ageo_ctrl = area_mean('control_qflux', 'fv_ageo', region='EAsia')
    easia_ageo_no_TIP = area_mean('no_TIP', 'fv_ageo', region='EAsia')
    easia_geo_ctrl = area_mean('control_qflux', 'dphidx', region='EAsia')
    easia_geo_no_TIP = area_mean('no_TIP', 'dphidx', region='EAsia')
    
    fig = plt.figure()
    ax = plt.subplot(111)
    easia_ageo_ctrl.plot(color='k', linewidth=2)
    (-1.*easia_geo_ctrl).plot(color='k', linestyle='--', linewidth=2)
    easia_ageo_no_TIP.plot(color='r', linewidth=2)
    (-1.*easia_geo_no_TIP).plot(color='r', linestyle='--', linewidth=2)
    plt.xlabel('Pentad')
    plt.ylabel('Momentum tendency, m/s/day')
    plt.xlim([1,72])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    #ax.legend(['Control', 'No Tibet'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
    plt.subplots_adjust(left=0.15, right=0.75, top=0.95, bottom=0.2)
    plt.savefig(plot_dir + 'fv_easia_no_TIP.pdf', format='pdf')
    plt.close()
    
    
    
    
    