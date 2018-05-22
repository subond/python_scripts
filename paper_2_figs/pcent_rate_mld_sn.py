"""
Produce plots showing precip centroid lat vs rate for experiments with varying year length, rotation rate, and mixed layer depth 27/03/2018

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


def p_cent_rate(data, days=False):
    """
    Inputs:
        data - xarray dataset climatology including precipitation as either precipitation or as convection_rain and condensation_rain
        days - instructs if data is in daily or pentad means
    Returns:
        dpcentdt - rate of movement of the precipitation centroid
        dpcentdt2 - rate of acceleration of the precipitation centroid
    """
    # Get total precip
    try:
        data['precipitation'] = data.condensation_rain + data.convection_rain
    except:
        data['precipitation'] = data.precipitation
    data['precipitation'] = make_sym(data.precipitation)
    
    # Locate precipitation centroid
    precip_centroid(data)
    
    # If units in are days rather than pentads (e.g. for shorter years) convert to pentads
    if days:
        dpcentdt = gr.ddt(data.p_cent, secperunit = 86400.) * 86400.
    else:
        dpcentdt = gr.ddt(data.p_cent) * 86400. 
    dpcentdt2 = gr.ddt(dpcentdt, secperunit = 86400.) * 86400.
    
    return dpcentdt, dpcentdt2


def p_cent_grad_scatter(run, days=False, period_fac=1., ax=None, color='k', linewidth=1.):
    """
    Inputs:
        run - name of a run stored in the climatologies folder
        days - to pass to p_cent_rate, instructs if data is in daily or pentad means
        period_fac - optional input which should be the ratio of the climatology year to a 360 day year
                     if set, the rate will be scaled by this factor
        ax - optionally specify an axis to plot to
        color - color to plot in (default black)
        linewidth - width of line to plot (default 1.)
    Outputs:
        plots to specified/current axes, this can be saved and modified outside the function
    """
    # Open dataset
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')

    dpcentdt, dpcentdt2 = p_cent_rate(data, days) # Get rate of movement of precip centroid
        
    if ax == None:
        plt.plot(data.p_cent, dpcentdt*period_fac, color, linewidth=linewidth)
        #plt.plot([data.p_cent[-1],data.p_cent[0]], [dpcentdt[-1]*period_fac,dpcentdt[0]*period_fac], color, linewidth=linewidth)
    else:
        ax.plot(data.p_cent, dpcentdt*period_fac, color, linewidth=linewidth) # Plot scatter of precip centroid lat vs rate
        #ax.plot([data.p_cent[-1],data.p_cent[0]], [dpcentdt[-1]*period_fac,dpcentdt[0]*period_fac], color, linewidth=linewidth)

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
    
    
def set_plot_features(ax, title='', legend_labels=[], fontsize=10, leg_title=None):
    """
    Inputs:
        ax - axis to modify
        title - title for axis (default empty string)
        legend_labels - label for legend (default empty list)
    """
    # Shrink current axis by 10%
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    legend = ax.legend(legend_labels, loc='upper left', borderaxespad=0., fontsize=fontsize, title=leg_title, ncol=2) #bbox_to_anchor=(1.05, 1),
    legend.get_title().set_fontsize(fontsize)
    ax.set_xlim([-25,25])
    ax.set_ylim([-1., 1.]) 
    #ax.set_title(title, fontsize=16)
    ax.grid(True,linestyle=':')



if __name__ == "__main__":
    
    # Set plotting directory
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)

    # Set figure parameters
    rcParams['figure.figsize'] = 15, 10.5
    rcParams['font.size'] = 14

    # Start figure with 4 subplots
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3)

    # Set colors for lines
    colors=['b','g','r','c','m','y']
    
    '''Mixed layer depths'''
    runs = ['mld_2.5', 'mld_5', 'sn_1.000', 
            'mld_15', 'mld_20']
    mlds = np.array([2.5,5.,10.,15.,20.])
    i=0
    max_rate_mld, max_rate_lat_mld, max_lat_mld = p_cent_rate_max(runs)
    dpdt_eq_mld = rate_at_eq(runs)
    for run in runs:
        if run == 'sn_1.000':
            p_cent_grad_scatter(run, color='k', ax=ax3, linewidth=2.)    
        else:
            p_cent_grad_scatter(run, color=colors[i], ax=ax3, linewidth=1.)
            i=i+1
    set_plot_features(ax3, title='Varying MLDs', legend_labels=['2.5', '5.', '10.', '15.', '20.'], leg_title='MLD (m)')
    
    A = np.array([ mlds, np.ones(mlds.shape) ])
    model = sm.OLS(max_rate_mld.values, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('=== Coeffs ===')
    print(consts[0], consts[1])
    print('=== Std Errs ===')
    print(2.*std_err[0], 2*std_err[1])
    line = mlds*consts[0] + consts[1]
    ax6.plot(mlds, max_rate_mld.values, 'xk', mew=2, ms=10)
    ax6.plot(mlds, line,'k')
    ax6.set_ylabel('ITCZ migration rate')
    
    model = sm.OLS(max_lat_mld.values, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('=== Coeffs ===')
    print(consts[0], consts[1])
    print('=== Std Errs ===')
    print(2.*std_err[0], 2*std_err[1])
    line = mlds*consts[0] + consts[1]
    ax9.plot(mlds, max_lat_mld.values, 'xk', mew=2, ms=10)
    ax9.plot(mlds, max_rate_lat_mld.values, 'xr', mew=2, ms=10)
    ax9.plot(mlds, line,'k')
    ax9.set_ylabel('Latitude')
    ax9.set_xlabel('MLD, m')
    
    A = np.array([ 1./mlds, np.ones(mlds.shape) ])
    model = sm.OLS(dpdt_eq_mld, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('=== Coeffs ===')
    print(consts[0], consts[1])
    print('=== Std Errs ===')
    print(2.*std_err[0], 2*std_err[1])
    line = 1./np.arange(2.5,20.01,0.01)*consts[0] + consts[1]
    ax6.plot(mlds, dpdt_eq_mld, 'xr', mew=2, ms=10)
    ax6.plot(np.arange(2.5,20.01,0.01), line,'r')
    ax6.set_xlabel('MLD, m')
    
    
    # Year length
    runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
            'sn_2.000', 'sn_3.000', 'sn_4.000']
    days = [True, False, False, False, False, False]
    max_rate_sn, max_rate_lat_sn, max_lat_sn = p_cent_rate_max(runs, days=days)
    dpdt_eq_sn = rate_at_eq(runs, days=days)
    # Set scaling factors for rate
    period_fac = [0.25, 0.5, 1., 2., 3., 4.]
    # Do 0.25 first as this is in days
    p_cent_grad_scatter(runs[0], days=True, color=colors[0], ax=ax1, linewidth=1.)
    # Do other plots, plot sn_1.000 run as thicker line
    j=1
    for i in range(1,6):
        if runs[i] == 'sn_1.000':
            p_cent_grad_scatter(runs[i], color='k', ax=ax1, linewidth=2.)    
        else:
            p_cent_grad_scatter(runs[i], color=colors[j], ax=ax1, linewidth=1.) 
            j=j+1
    set_plot_features(ax1, title='Varying Year Lengths', legend_labels=['0.25', '0.5', '1.0', '2.0', '3.0', '4.0'], leg_title='P/P$_{E}$')
    
    
    ax4.plot(period_fac, max_rate_sn.values, 'xk', mew=2, ms=10)
    ax4.plot(period_fac, dpdt_eq_sn, 'xr', mew=2, ms=10)
    ax4.set_ylabel('ITCZ migration rate')
    ax4.set_xlabel('P/P$_{E}$')
    #ax5.set_xscale('log')
    #ax5.set_yscale('log')
    
    ax7.plot(period_fac, max_lat_sn.values, 'xk', mew=2, ms=10)
    ax7.plot(period_fac, max_rate_lat_sn.values, 'xr', mew=2, ms=10)
    ax7.set_ylabel('Latitude')
    ax7.set_xlabel('P/P$_{E}$')
    #ax6.set_xscale('log')
    #ax6.set_yscale('log')
    
    
    # Year length, scaled
    runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
            'sn_2.000', 'sn_3.000', 'sn_4.000']
    # Set scaling factors for rate
    period_fac = [0.25, 0.5, 1., 2., 3., 4.]
    # Do 0.25 first as this is in days
    p_cent_grad_scatter(runs[0], days=True, period_fac=period_fac[0], color=colors[0], ax=ax2, linewidth=1.)
    # Do other plots, plot sn_1.000 run as thicker line
    j=1
    for i in range(1,6):
        if runs[i] == 'sn_1.000':
            p_cent_grad_scatter(runs[i], color='k', ax=ax2, linewidth=2.)    
        else:
            p_cent_grad_scatter(runs[i], period_fac=period_fac[i], color=colors[j], ax=ax2, linewidth=1.)    
            j=j+1
    set_plot_features(ax2, title='Varying Year Lengths (Rates scaled)', legend_labels=['0.25', '0.5', '1.0', '2.0', '3.0', '4.0'], leg_title='P/P$_{E}$')
    
    ax5.plot(period_fac, max_rate_sn.values[:,0]*period_fac, 'xk', mew=2, ms=10)
    ax5.plot(period_fac, dpdt_eq_sn*period_fac, 'xr', mew=2, ms=10)
    ax5.set_xlabel('P/P$_{E}$')
    ax5.set_ylabel('ITCZ migration rate')
    
    
    
    for ax in [ax1, ax2, ax3]:
        ax.set_ylabel('P. cent. lat. rate')
        ax.set_xlabel('Precip centroid lat.')

    for ax in [ax4, ax5, ax6]:
        ax.set_ylim([0,1])
    
    for ax in [ax7, ax8]:
        ax.set_ylim([0,25])
        
    ax8.axis('off') 
    
    plt.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.05, hspace=0.3, wspace=0.3)
    
    # Save as a pdf
    plt.savefig(plot_dir + 'pcent_rate_mld_sn.pdf', format='pdf')
    plt.close()
    
    