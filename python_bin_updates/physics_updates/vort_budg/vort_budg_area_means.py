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
from plot_vorticity_breakdown import vort_budg_terms


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
    
    data = vort_budg_terms(run, ll=True)
    
    coslat = np.cos(data[var].lat * np.pi/180.)
    sinlat = np.sin(data[var].lat * np.pi/180.)
    
    var_area_weighted = data[var].sel(pfull=150) * coslat
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= latin[0] and data.lat[i] < latin[1]]
    lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    
    var_mean = var_area_weighted.sel(lon=lons).sel(lat=lats).sum(('lat','lon')) / (coslat.sel(lat=lats).sum('lat') * len(lons))

    return var_mean


        

if __name__ == "__main__":
    
    rcParams['figure.figsize'] = 5, 3
    rcParams['font.size'] = 14
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/climatology/mom_budg_area_means/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #india_vor_ctrl = area_mean('control_qflux', 'vor', region='India')
    #india_vor_no_TIP = area_mean('no_TIP', 'vor', region='India')
    
    #fig = plt.figure()
    #ax = plt.subplot(111)
    #india_vor_ctrl.plot(color='k', linewidth=2)
    #india_vor_no_TIP.plot(color='r', linewidth=2)
    #plt.xlabel('Pentad')
    #plt.ylabel('Absolute vorticity, /day')
    #plt.xlim([1,72])
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    #ax.legend(['Control', 'No Tibet'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
    #plt.subplots_adjust(left=0.15, right=0.75, top=0.95, bottom=0.2)
    #plt.savefig(plot_dir + 'vor_india_no_TIP.pdf', format='pdf')
    #plt.close()
    
    
    easia_vor_ctrl = area_mean('control_qflux', 'vor', region='EAsia')
    easia_vor_no_TIP = area_mean('no_TIP', 'vor', region='EAsia')
    easia_vor_no_am = area_mean('no_americas', 'vor', region='EAsia')
    
    fig = plt.figure()
    ax = plt.subplot(111)
    easia_vor_ctrl.plot(color='k', linewidth=2)
    #easia_vor_no_TIP.plot(color='r', linewidth=2)
    plt.title('')
    plt.xlabel('Pentad')
    plt.ylabel('Absolute vorticity, day$^{-1}$')
    plt.xlim([1,72])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    #ax.legend(['Control', 'No Tibet'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
    plt.subplots_adjust(left=0.15, right=0.75, top=0.95, bottom=0.2)
    plt.savefig(plot_dir + 'vor_easia.pdf', format='pdf')
    plt.close()
    
    fig = plt.figure()
    ax = plt.subplot(111)
    easia_vor_ctrl.plot(color='k', linewidth=2)
    easia_vor_no_TIP.plot(color='r', linewidth=2)
    plt.title('')
    plt.xlabel('Pentad')
    plt.ylabel('Absolute vorticity, day$^{-1}$')
    plt.xlim([1,72])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    #ax.legend(['Control', 'No Tibet'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
    plt.subplots_adjust(left=0.15, right=0.75, top=0.95, bottom=0.2)
    plt.savefig(plot_dir + 'vor_easia_no_TIP.pdf', format='pdf')
    plt.close()
    
    
    fig = plt.figure()
    ax = plt.subplot(111)
    easia_vor_ctrl.plot(color='k', linewidth=2)
    easia_vor_no_am.plot(color='r', linewidth=2)
    plt.title('')
    plt.xlabel('Pentad')
    plt.ylabel('Absolute vorticity, day$^{-1}$')
    plt.xlim([1,72])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    #ax.legend(['Control', 'No Tibet'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
    plt.subplots_adjust(left=0.15, right=0.75, top=0.95, bottom=0.2)
    plt.savefig(plot_dir + 'vor_easia_no_am.pdf', format='pdf')
    plt.close()
    
    
    
    
    