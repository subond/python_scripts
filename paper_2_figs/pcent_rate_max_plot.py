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
from pcent_rate_max import p_cent_rate_max
import statsmodels.api as sm


if __name__ == "__main__":
    
    # Set plotting directory
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    

    rcParams['figure.figsize'] = 5, 6
    rcParams['font.size'] = 16
    
    runs = ['rt_0.500', 'rt_0.750', 'sn_1.000',
            'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']
    
    #runs_5 = ['rt_0.500_5','rt_0.750_5', 'mld_5', 
    #         'rt_1.250_5','rt_1.500_5','rt_1.750_5','rt_2.000_5']
    
    #runs_15 = ['rt_0.500_15','rt_0.750_15', 'mld_15', 
    #         'rt_1.250_15','rt_1.500_15','rt_1.750_15','rt_2.000_15']
    
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    max_rate, max_rate_lat, max_lat = p_cent_rate_max(runs)
    #max_rate_5,  max_rate_lat_5, max_lat_5 = p_cent_rate_max(runs_5)
    #max_rate_15, max_rate_lat_15, max_lat_15 = p_cent_rate_max(runs_15)
    
    rots = np.array([0.5,0.75,1.,1.25,1.5,1.75,2.])
    A = np.array([ np.log(rots), np.ones(rots.shape) ])
    model = sm.OLS(np.log(max_rate_lat.values), A.T)
    #A = np.array([ rots, np.ones(rots.shape) ])
    #model = sm.OLS(max_rate_lat.values, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('=== Coeffs ===')
    print(np.exp(consts[1]), consts[0])
    #print(consts[0], consts[1])
    print('=== Std Errs ===')
    print(2.*std_err[1]*np.exp(consts[1]), 2*std_err[0])
    #print(2.*std_err[0], 2*std_err[1])
    
    line = np.exp(consts[1]) * rots**(consts[0])
    #line = consts[1] + rots * consts[0]
    
        
    ax1.plot([0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.], max_rate_lat.values, 'xk', mew=2, ms=10)
    #ax1.plot([0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.], max_rate_lat_5.values, 'xb', mew=2, ms=10)
    #ax1.plot([0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.], max_rate_lat_15.values, 'xr', mew=2, ms=10)
    ax1.plot(rots, line,'k')
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')
    ax1.set_ylabel('Lat. of max rate')
    
    A = np.array([ np.log(rots), np.ones(rots.shape) ])
    #A = np.array([ rots**(-2./3.), np.ones(rots.shape) ])
    model = sm.OLS(np.log(max_lat.values), A.T)
    #model = sm.OLS(max_lat.values, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('=== Coeffs ===')
    print(np.exp(consts[1]), consts[0])
    print('=== Std Errs ===')
    print(2*std_err[1]*np.exp(consts[1]), 2*std_err[0])
    #print(2.*std_err[0], 2*std_err[1])
    
    line = np.exp(consts[1]) * rots**(consts[0])
    #line = consts[1] + consts[0]*rots**(-2./3.)
    
    ax2.plot([0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.], max_lat.values, 'xk', mew=2, ms=10)
    ax2.plot(rots, line,'k')
    #ax2.plot([0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.], max_lat_5.values, 'xb', mew=2, ms=10)
    #ax2.plot([0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.], max_lat_15.values, 'xr', mew=2, ms=10)
    ax2.set_xlabel('')
    #ax2.set_ylim([0,1.1])
    #ax2.set_yticks([0,0.25,0.5,0.75])
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')
    ax2.set_ylabel('Amplitude')
    
    ax2.set_xlabel('$\Omega$/$\Omega_{E}$')
    
    plt.subplots_adjust(right=0.95, left=0.2, top=0.95, bottom=0.1, hspace=0.15)
    
    plt.savefig(plot_dir + 'rotation_scatter.pdf', format='pdf')
    plt.close()
    
    
    
    