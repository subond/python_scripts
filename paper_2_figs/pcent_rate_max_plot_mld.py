"""
Calculate rate of movement of precipitation centroid for different seasons 02/02/2018

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_centroid
from data_handling_updates import gradients as gr, make_sym
from pylab import rcParams
from pcent_rate_max import p_cent_rate_max
import statsmodels.api as sm

def rate_at_eq(runs, do_make_sym=True, days=None):
    dpdt_eq = []
    if days==None:
        days=[False]*len(runs)
    j=0
    for run in runs: 
        data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')

        if do_make_sym: # If symmetric data is wanted, average over both hemispheres (NB currently only set up for a climatology 2/05/18)
            data['precipitation'] = make_sym(data.precipitation)

        # Locate precipitation centroid
        precip_centroid(data)
        
        # Get rate of movement of precip centroid
        if days[j]:
            dpcentdt = gr.ddt(data.p_cent, secperunit = 86400.) * 86400.
        else:
            dpcentdt = gr.ddt(data.p_cent) * 86400.
        #dpcentdt_max = dpcentdt_pmask.where(dpcentdt_pmask==dpcentdt_pmask.max('xofyear'),drop=True)   # Find the maximum rate
        p_cent = np.abs(data.p_cent.where(dpcentdt>=0.))        
        dpdt_eq_j = dpcentdt.where(p_cent == p_cent.min('xofyear'), drop=True)
        dpdt_eq.append(dpdt_eq_j.values[0])
        j=j+1
    return np.asarray(dpdt_eq)
    

if __name__ == "__main__":
    
    # Set plotting directory
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    

    rcParams['figure.figsize'] = 5, 6
    rcParams['font.size'] = 16
    
    runs = ['mld_2.5', 'mld_5', 'sn_1.000', 'mld_15', 'mld_20']
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    max_rate, max_rate_lat, max_lat = p_cent_rate_max(runs)
    dpdt_eq = rate_at_eq(runs)
    
    mlds = np.array([2.5,5.,10.,15.,20.])
    A = np.array([ mlds, np.ones(mlds.shape) ])
    
    model = sm.OLS(max_rate.values, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('=== Coeffs ===')
    print(consts[0], consts[1])
    print('=== Std Errs ===')
    print(2.*std_err[0], 2*std_err[1])
    line = mlds*consts[0] + consts[1]
    
    ax1.plot(mlds, max_rate.values, 'xk', mew=2, ms=10)
    ax1.plot(mlds, line,'k')
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')
    ax1.set_ylabel('Max rate')
    
    
    model = sm.OLS(max_lat.values, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('=== Coeffs ===')
    print(consts[0], consts[1])
    print('=== Std Errs ===')
    print(2.*std_err[0], 2*std_err[1])
    line = mlds*consts[0] + consts[1]
    
    ax2.plot(mlds, max_lat.values, 'xk', mew=2, ms=10)
    ax2.plot(mlds, line,'k')
    #ax2.set_ylim([0,1.1])
    #ax2.set_yticks([0,0.25,0.5,0.75])
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')
    ax2.set_ylabel('Amplitude')
    
    A = np.array([ 1./mlds, np.ones(mlds.shape) ])
    model = sm.OLS(dpdt_eq, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('=== Coeffs ===')
    print(consts[0], consts[1])
    print('=== Std Errs ===')
    print(2.*std_err[0], 2*std_err[1])
    line = 1./mlds*consts[0] + consts[1]
    
    ax3.plot(mlds, dpdt_eq, 'xk', mew=2, ms=10)
    ax3.plot(mlds, line,'k')
    #ax2.set_ylim([0,1.1])
    #ax2.set_yticks([0,0.25,0.5,0.75])
    #ax3.set_xscale('log')
    #ax3.set_yscale('log')
    ax3.set_ylabel('Rate at Eq.')
    
    ax3.set_xlabel('MLD, m')
    
    plt.subplots_adjust(right=0.95, left=0.2, top=0.95, bottom=0.1, hspace=0.15)
    
    plt.savefig(plot_dir + 'mld_scatter.pdf', format='pdf')
    plt.close()
    
    
    
    
    