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
    rcParams['figure.figsize'] = 10, 9
    rcParams['font.size'] = 12
    
    subsolar_point = -23.439*np.cos(np.arange(0.,361.,1.)*np.pi/180.)
    subsolar_point = xr.DataArray(subsolar_point, coords=[np.arange(0.,361.,1.)], dims=['xofyear'])
    subsolar_rate = gr.ddt(subsolar_point, secperunit=86400.) * 86400.    
    
    # Start figure with 4 subplots
    #fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3)
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
    
    # Switch on errorbars
    errorbars=True
    
    # Set colors for lines
    #colors=['b','g','r','c','m','y']
    colors=['C0','C1','C2','C3','C4','C5']
    
    
    '''Seasons'''
    runs_sn = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
            'sn_2.000', 'sn_3.000', 'sn_4.000']
    days = [True, False, False, False, False, False]
    period_fac = [0.25, 0.5, 1., 2., 3., 4.]
    
    # Get rates and lats needed
    max_rate_sn, max_rate_lat_sn, max_lat_sn = p_cent_rate_max(runs_sn, days=days)
    dpdt_eq_sn = rate_at_eq(runs_sn, days=days)
    
    # Print values needed in paper
    print('sn1 max rate: ', max_rate_sn[2].values)
    print('sn1 max_rate lat: ', max_rate_lat_sn[2].values)
    print('sn1 eq rate: ', dpdt_eq_sn[2])
    print('scaled rates at eq: ', dpdt_eq_sn * period_fac)
    
    # Calculate 95% confidence interval from bootstrapping data
    err_sn=[]
    for errname in runs_sn:
        err_sn.append(np.load(errname+'_bootstrap.npy'))
    err_sn = np.asarray(err_sn)
    lower_sn = np.percentile(err_sn,2.5,axis=2)
    upper_sn = np.percentile(err_sn,97.5,axis=2)
    
    # Plot 0.25 first as this is in days
    p_cent_grad_scatter(runs_sn[0], days=True, color=colors[0], ax=ax1, linewidth=1.)
    # Do other plots, plot sn_1.000 run as thicker line
    j=1
    for run in runs_sn[1:]:
        if run == 'sn_1.000':
            p_cent_grad_scatter(run, color='k', ax=ax1, linewidth=2.)    
        else:
            p_cent_grad_scatter(run, color=colors[j], ax=ax1, linewidth=1.) 
            j=j+1
    set_plot_features(ax1, title='Varying Year Lengths', legend_labels=['0.25', '0.5', '1.0', '2.0', '3.0', '4.0'], leg_title='P/P$_{E}$')
    
    
    # Now do unscaled rate plots and lat plots
    
    if errorbars:
        ax3.errorbar(period_fac, max_rate_sn.values,
                 yerr=[max_rate_sn.values-lower_sn[:,0], upper_sn[:,0]-max_rate_sn.values], 
                 linestyle='none', color='k', marker='.',mew=2, ms=8)
    
        ax3.errorbar(period_fac, dpdt_eq_sn,
                 yerr=[dpdt_eq_sn-lower_sn[:,3], upper_sn[:,3]-dpdt_eq_sn], 
                 linestyle='none', color='r', marker='.',mew=2, ms=8)
    
        ax5.errorbar(period_fac, max_lat_sn.values,
                 yerr=[max_lat_sn.values-lower_sn[:,2], upper_sn[:,2]-max_lat_sn.values], 
                 linestyle='none', color='r', marker='.',mew=2, ms=8)
    
        ax5.errorbar(period_fac, max_rate_lat_sn.values,
                 yerr=[max_rate_lat_sn.values-lower_sn[:,1], upper_sn[:,1]-max_rate_lat_sn.values], 
                 linestyle='none', color='k', marker='.',mew=2, ms=8)
    
    else:
        ax3.plot(period_fac, max_rate_sn.values, 'xk', mew=2, ms=10)
        ax3.plot(period_fac, dpdt_eq_sn, 'xr', mew=2, ms=10)
        ax5.plot(period_fac, max_lat_sn.values, 'xr', mew=2, ms=10)
        ax5.plot(period_fac, max_rate_lat_sn.values, 'xk', mew=2, ms=10)
    
    # Set labels
    ax3.set_ylabel('ITCZ migration rate')
    ax3.set_xlabel('P/P$_{E}$')
    ax5.set_ylabel('Latitude')
    ax5.set_xlabel('P/P$_{E}$')

    
    
    '''Seasons scaled'''

    # Plot 0.25 first as this is in days
    p_cent_grad_scatter(runs_sn[0], days=True, period_fac=period_fac[0], color=colors[0], ax=ax2, linewidth=1.)
    # Do other plots, plot sn_1.000 run as thicker line
    j=1
    period_fac_no1 = [0.25, 0.5, 2., 3., 4.]
    for run in runs_sn[1:]:
        if run == 'sn_1.000':
            p_cent_grad_scatter(run, color='k', ax=ax2, linewidth=2.)    
        else:
            p_cent_grad_scatter(run, period_fac=period_fac_no1[j], color=colors[j], ax=ax2, linewidth=1.)    
            j=j+1
    set_plot_features(ax2, title='Varying Year Lengths (Rates scaled)', legend_labels=['0.25', '0.5', '1.0', '2.0', '3.0', '4.0'], leg_title='P/P$_{E}$')
    
    # Now do scaled rate plots and lat plots
    
    if errorbars:
        ax4.errorbar(period_fac, max_rate_sn.values[:,0]*period_fac,
                 yerr=[(max_rate_sn.values[:,0]-lower_sn[:,0])*period_fac, (upper_sn[:,0]-max_rate_sn.values[:,0])*period_fac], 
                 linestyle='none', color='k', marker='.',mew=2, ms=8)
    
        ax4.errorbar(period_fac, dpdt_eq_sn*period_fac,
                 yerr=[(dpdt_eq_sn-lower_sn[:,3])*period_fac, (upper_sn[:,3]-dpdt_eq_sn)*period_fac], 
                 linestyle='none', color='r', marker='.',mew=2, ms=8)
    
    else:
        ax4.plot(period_fac, max_rate_sn.values[:,0]*period_fac, 'xk', mew=2, ms=10)
        ax4.plot(period_fac, dpdt_eq_sn*period_fac, 'xr', mew=2, ms=10)
    
    # Set labels
    ax4.set_xlabel('P/P$_{E}$')
    ax4.set_ylabel('Scaled ITCZ migration rate')
    
    ax1.text(-35, 1., 'a)')
    ax2.text(-37, 1., 'd)')
    ax3.text(-0.7, 1., 'b)')
    ax4.text(-0.9, 1., 'e)')
    ax5.text(-0.7, 25., 'c)')
    
    ax1.set_ylabel('ITCZ migration rate')
    ax2.set_ylabel('Scaled ITCZ migration rate')
    
    
    for ax in [ax1, ax2]:
        ax.set_xlabel('ITCZ latitude')

    for ax in [ax3, ax4]:
        ax.set_ylim([0,1])
        ax.grid(True,linestyle=':')
    
    for ax in [ax5]:
        ax.set_ylim([0,25])
        ax.grid(True,linestyle=':')
        
    ax6.axis('off') 
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.96, bottom=0.07, hspace=0.3, wspace=0.3)
    
    # Save as a pdf
    plt.savefig(plot_dir + 'pcent_rate_sn.pdf', format='pdf')
    plt.close()
    
    
    
    
    
    
    '''Mixed layer depths'''
    
    # Set figure parameters
    rcParams['figure.figsize'] = 5, 9
    rcParams['font.size'] = 12

    fig, ((ax1, ax2, ax3)) = plt.subplots(3, 1)
    
    
    runs_mld = ['mld_2.5', 'mld_5', 'sn_1.000', 
            'mld_15', 'mld_20']
    mlds = np.array([2.5,5.,10.,15.,20.])
    
    # Get rates and lats needed
    max_rate_mld, max_rate_lat_mld, max_lat_mld = p_cent_rate_max(runs_mld)
    dpdt_eq_mld = rate_at_eq(runs_mld)
    
    # Do plots, plot sn_1.000 run as thicker line
    j=0
    for run in runs_mld:
        if run == 'sn_1.000':
            p_cent_grad_scatter(run, color='k', ax=ax1, linewidth=2.)    
        else:
            p_cent_grad_scatter(run, color=colors[j], ax=ax1, linewidth=1.)
            j=j+1
    set_plot_features(ax1, title='Varying MLDs', legend_labels=['2.5', '5.', '10.', '15.', '20.'], leg_title='Mixed layer depth, m')
    
    #ax1.plot(subsolar_point, subsolar_rate, 'm', linewidth=2.)
    
    print('mld peak lat:', np.mean(max_rate_lat_mld[0:4].values), '+-', np.std(max_rate_lat_mld[0:4].values))
    print('sn peak lat:', np.mean(max_rate_lat_sn[1:6].values), '+-', np.std(max_rate_lat_sn[1:6].values))
    
    # Calculate 95% confidence interval from bootstrapping data
    err_mld=[]
    for errname in runs_mld:
        err_mld.append(np.load(errname+'_bootstrap.npy'))
    err_mld = np.asarray(err_mld)
    lower_mld = np.percentile(err_mld,2.5,axis=2)
    upper_mld = np.percentile(err_mld,97.5,axis=2)
    
    if errorbars:
        ax2.errorbar(mlds, max_rate_mld.values,
                     yerr=[max_rate_mld.values-lower_mld[:,0], upper_mld[:,0]-max_rate_mld.values], 
                     linestyle='none', color='k', marker='.',mew=2, ms=8)
    
        ax2.errorbar(mlds, dpdt_eq_mld,
                     yerr=[dpdt_eq_mld-lower_mld[:,3], upper_mld[:,3]-dpdt_eq_mld], 
                     linestyle='none', color='r', marker='.',mew=2, ms=8)
        
        ax3.errorbar(mlds, max_lat_mld.values,
                     yerr=[max_lat_mld.values-lower_mld[:,2], upper_mld[:,2]-max_lat_mld.values], 
                     linestyle='none', color='r', marker='.',mew=2, ms=8)
    
        ax3.errorbar(mlds, max_rate_lat_mld.values,
                     yerr=[max_rate_lat_mld.values-lower_mld[:,1], upper_mld[:,1]-max_rate_lat_mld.values], 
                     linestyle='none', color='k', marker='.',mew=2, ms=8)
    
    else:
        ax2.plot(mlds, max_rate_mld.values, 'xk', mew=2, ms=10)
        ax2.plot(mlds, dpdt_eq_mld, 'xr', mew=2, ms=10)
        ax3.plot(mlds, max_lat_mld.values, 'xr', mew=2, ms=10)
        ax3.plot(mlds, max_rate_lat_mld.values, 'xk', mew=2, ms=10)
        
    # Fit a linear relation to the max rates
    A = np.array([ mlds, np.ones(mlds.shape) ])
    model = sm.OLS(max_rate_mld.values, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('* MLD max rate = (', consts[0], ' +- ', 2*std_err[0], ') * MLD + (', consts[1], ' +- ', 2*std_err[1], ') *')
    line_maxrate = mlds*consts[0] + consts[1]
    
    
    # Fit a linear relation to the max lats
    model = sm.OLS(max_lat_mld.values, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('* MLD max lat = (', consts[0], ' +- ', 2*std_err[0], ') * MLD + (', consts[1], ' +- ', 2*std_err[1], ') *')
    line_maxlat = mlds*consts[0] + consts[1]
    
    # Fit an inverse relation to the rates at the equator
    A = np.array([ 1./mlds, np.ones(mlds.shape) ])
    model = sm.OLS(dpdt_eq_mld, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    print('* MLD eq rate = (', consts[0], ' +- ', 2*std_err[0], ') / MLD + (', consts[1], ' +- ', 2*std_err[1], ') *')
    line_eqrate = 1./np.arange(2.5,20.01,0.01)*consts[0] + consts[1]
    
    
    # Plot the best fit line and add labels        
    ax2.plot(mlds, line_maxrate,'k')
    ax2.plot(np.arange(2.5,20.01,0.01), line_eqrate,'r')
    ax2.set_ylabel('ITCZ migration rate')
    ax2.set_xlabel('Mixed layer depth, m')
    
    ax3.plot(mlds, line_maxlat,'r')
    ax3.set_ylabel('Latitude')
    ax3.set_xlabel('Mixed layer depth, m')
        
    
    ax1.text(-35, 1., 'a)')
    ax2.text(-2.5, 1., 'b)')
    ax3.text(-2.5, 25., 'c)')
    
    ax1.set_ylabel('ITCZ migration rate')
    ax1.set_xlabel('ITCZ latitude')

    ax2.set_ylim([0,1])
    ax2.grid(True,linestyle=':')
    
    ax3.set_ylim([0,25])
    ax3.grid(True,linestyle=':')
        
    
    plt.subplots_adjust(left=0.2, right=0.97, top=0.96, bottom=0.07, hspace=0.3, wspace=0.3)
    
    # Save as a pdf
    plt.savefig(plot_dir + 'pcent_rate_mld.pdf', format='pdf')
    plt.close()
    
    
    
    