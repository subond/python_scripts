"""
Calculate rate of movement of precipitation centroid for different seasons 02/02/2018

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_centroid
from data_handling_updates import gradients as gr
from pylab import rcParams
from data_handling_updates import make_sym


def p_cent_rate_max(runs, days=None):
    # Get the maximum rate of change of the precipitation centroid latitude, and the latitude at which this occurs.
    
    max_rate = []
    max_rate_lat = []
    if days==None:
        days=[False]*len(runs)
        
    j=0
    for run in runs:
        # Open dataset
        data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
        
        # Get total precip
        try:
            data['precipitation'] = data.condensation_rain + data.convection_rain
        except:
            data['precipitation'] = data.precipitation
        
        data['precipitation'] = make_sym(data.precipitation)

        # Locate precipitation centroid
        precip_centroid(data)
            
        # Get rate of movement of precip centroid
        if days[j]:
            dpcentdt = gr.ddt(data.p_cent, secperunit = 86400.) * 86400.
        else:
            dpcentdt = gr.ddt(data.p_cent) * 86400.
                    
    
        dpcentdt_ppos = dpcentdt.where(data.p_cent>=0.)   # Find precip centroid rate where precip centroid is in the northern hemisphere
        dpcentdt_max = dpcentdt_ppos.where(dpcentdt_ppos==dpcentdt_ppos.max(),drop=True)   # Find the maximum of the above
        if len(dpcentdt_max) > 1:
            dpcentdt_max = dpcentdt_max[0]
        pcent_dtmax = data.p_cent.sel(xofyear=dpcentdt_max.xofyear)    # Find the location of the preciptiation when the rate is maximum
        
        print(dpcentdt_max.values, pcent_dtmax.values)     # Print both
        
        max_rate.append(dpcentdt_max)
        max_rate_lat.append(pcent_dtmax)
        
        j=j+1
    
    max_rate = np.asarray(max_rate)
    max_rate_lat = np.asarray(max_rate_lat)
    
    return max_rate, max_rate_lat
    
    

if __name__ == "__main__":
    
    # Set plotting directory
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    

    rcParams['figure.figsize'] = 5, 6
    rcParams['font.size'] = 16
    
    runs = ['rt_0.500', 'rt_0.750', 'sn_1.000',
            'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']
    
    runs_5 = ['rt_0.750_5', 'mld_5', 'rt_1.250_5']
    
    runs_15 = ['rt_0.750_15', 'mld_15', 'rt_1.250_15']
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    max_rate, max_rate_lat = p_cent_rate_max(runs)
    max_rate_5,  max_rate_lat_5 = p_cent_rate_max(runs_5)
    max_rate_15, max_rate_lat_15 = p_cent_rate_max(runs_15)
    
    ax1.plot([0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.], max_rate_lat, 'xk', mew=2, ms=10)
    ax1.plot([0.75, 1., 1.25], max_rate_lat_5, 'xb', mew=2, ms=10)
    ax1.plot([0.75, 1., 1.25], max_rate_lat_15, 'xr', mew=2, ms=10)
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')
    ax1.set_ylabel('Lat. of max rate')
    
    ax2.plot([0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.], max_rate, 'xk', mew=2, ms=10)
    ax2.plot([0.75, 1., 1.25], max_rate_5, 'xb', mew=2, ms=10)
    ax2.plot([0.75, 1., 1.25], max_rate_15, 'xr', mew=2, ms=10)
    ax2.set_xlabel('')
    #ax2.set_ylim([0,1.1])
    ax2.set_yticks([0,0.25,0.5,0.75])
    ax2.set_ylabel('Max rate')
    
    ax2.set_xlabel('$\Omega$/$\Omega_{E}$')
    
    plt.subplots_adjust(right=0.95, left=0.2, top=0.95, bottom=0.1, hspace=0.1)
    
    plt.savefig(plot_dir + 'rotation_scatter.pdf', format='pdf')
    plt.close()
    
    
    
    
    