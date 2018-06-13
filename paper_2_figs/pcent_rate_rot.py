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
    
    
def set_plot_features(ax, title='', do_legend=True, legend_labels=[], fontsize=10, leg_title=None):
    """
    Inputs:
        ax - axis to modify
        title - title for axis (default empty string)
        legend_labels - label for legend (default empty list)
    """
    # Shrink current axis by 10%
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    if do_legend:
        legend = ax.legend(legend_labels, loc='upper left', borderaxespad=0., fontsize=fontsize, title=leg_title, ncol=2) #bbox_to_anchor=(1.05, 1),
        legend.get_title().set_fontsize(fontsize)
    ax.set_xlim([-30,30])
    ax.set_ylim([-1., 1.]) 
    ax.set_xticks(np.arange(-30,31,10)) 
    #ax.set_title(title, fontsize=16)
    ax.grid(True,linestyle=':')



if __name__ == "__main__":
    
    # Set plotting directory
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)

    # Set figure parameters
    rcParams['figure.figsize'] = 15, 7
    rcParams['font.size'] = 14
    
    errorbars=True
    
    
    ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
    ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
    ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
    ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
    ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)
    
    # Set colors for lines
    colors=['b','g','r','c','m','y']
    
    '''Rotation'''
    runs = ['rt_0.500', 'rt_0.750', 'sn_1.000',
            'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']
            
    runs_5 = ['rt_0.500_5', 'rt_0.750_5', 'mld_5',
            'rt_1.250_5', 'rt_1.500_5', 'rt_1.750_5', 'rt_2.000_5']
            
    runs_15 = ['rt_0.500_15', 'rt_0.750_15', 'mld_15',
            'rt_1.250_15', 'rt_1.500_15', 'rt_1.750_15', 'rt_2.000_15']
            
    rots = np.array([0.5,0.75,1.,1.25,1.5,1.75,2.])
    
    max_rate_rot, max_rate_lat_rot, max_lat_rot = p_cent_rate_max(runs)
    max_rate_rot_5, max_rate_lat_rot_5, max_lat_rot_5 = p_cent_rate_max(runs_5)
    max_rate_rot_15, max_rate_lat_rot_15, max_lat_rot_15 = p_cent_rate_max(runs_15)
    print('rot peak lat:', np.mean(max_rate_lat_rot.values), '+-', np.std(max_rate_lat_rot.values))
    
    i=0
    for run in runs_5:
        if run == 'mld_5':
            p_cent_grad_scatter(run, color='k', ax=ax1, linewidth=2.)    
        else:
            p_cent_grad_scatter(run, color=colors[i], ax=ax1, linewidth=1.)
            i=i+1
    set_plot_features(ax1, title='', do_legend=False) #legend_labels=['0.5', '0.75', '1.0', '1.25', '1.5', '1.75', '2.'], leg_title='$\Omega$/$\Omega_{E}$')
    ax1.set_ylabel('ITCZ migration rate')
    ax1.set_xlabel('Precip centroid lat.')
    
    i=0
    for run in runs:
        if run == 'sn_1.000':
            p_cent_grad_scatter(run, color='k', ax=ax2, linewidth=2.)    
        else:
            p_cent_grad_scatter(run, color=colors[i], ax=ax2, linewidth=1.)
            i=i+1
    set_plot_features(ax2, title='', do_legend=False) #, legend_labels=['0.5', '0.75', '1.0', '1.25', '1.5', '1.75', '2.'], leg_title='$\Omega$/$\Omega_{E}$')
    ax2.set_ylabel('ITCZ migration rate')
    ax2.set_xlabel('Precip centroid lat.')
    
    i=0
    for run in runs_15:
        if run == 'mld_15':
            p_cent_grad_scatter(run, color='k', ax=ax3, linewidth=2.)    
        else:
            p_cent_grad_scatter(run, color=colors[i], ax=ax3, linewidth=1.)
            i=i+1
    set_plot_features(ax3, title='', legend_labels=['0.5', '0.75', '1.0', '1.25', '1.5', '1.75', '2.'], leg_title='$\Omega$/$\Omega_{E}$')
    ax3.set_ylabel('ITCZ migration rate')
    ax3.set_xlabel('Precip centroid lat.')
    
    
    if errorbars:
        # Calculate 95% confidence interval from bootstrapping data and add to plots
        err=[]
        for errname in runs:
            err.append(np.load(errname+'_bootstrap.npy'))
        err = np.asarray(err)
        lower = np.percentile(err,2.5,axis=2)
        upper = np.percentile(err,97.5,axis=2)
        
        err_5=[]
        for errname in runs_5:
            err_5.append(np.load(errname+'_bootstrap.npy'))
        err_5 = np.asarray(err_5)
        lower_5 = np.percentile(err_5,2.5,axis=2)
        upper_5 = np.percentile(err_5,97.5,axis=2)
    
        err_15=[]
        for errname in runs_15:
            err_15.append(np.load(errname+'_bootstrap.npy'))
        err_15 = np.asarray(err_15)
        lower_15 = np.percentile(err_15,2.5,axis=2)
        upper_15 = np.percentile(err_15,97.5,axis=2)
        
        ax4.errorbar(rots, max_rate_rot_5.values,
                 yerr=[max_rate_rot_5.values-lower_5[:,0], upper_5[:,0]-max_rate_rot_5.values], 
                 linestyle='none', color='r', marker='.',mew=2, ms=8)

        ax4.errorbar(rots, max_rate_rot_15.values,
                 yerr=[max_rate_rot_15.values-lower_15[:,0], upper_15[:,0]-max_rate_rot_15.values], 
                 linestyle='none', color='b', marker='.',mew=2, ms=8)
    
        ax4.errorbar(rots, max_rate_rot.values,
                 yerr=[max_rate_rot.values-lower[:,0], upper[:,0]-max_rate_rot.values], 
                 linestyle='none', color='k', marker='.',mew=2, ms=8)
        
        ax5.errorbar(rots, max_rate_lat_rot_5.values,
                     yerr=[max_rate_lat_rot_5.values-lower_5[:,1], upper_5[:,1]-max_rate_lat_rot_5.values], 
                     linestyle='none', color='r', marker='.',mew=2, ms=8, label='5m')
                 
        ax5.errorbar(rots, max_rate_lat_rot_15.values,
                     yerr=[max_rate_lat_rot_15.values-lower_15[:,1], upper_15[:,1]-max_rate_lat_rot_15.values], 
                     linestyle='none', color='b', marker='.',mew=2, ms=8, label='15m')
    
        ax5.errorbar(rots, max_rate_lat_rot.values,
                     yerr=[max_rate_lat_rot.values-lower[:,1], upper[:,1]-max_rate_lat_rot.values], 
                     linestyle='none', color='k', marker='.',mew=2, ms=8, label='10m')
    
    else:
        # Otherwise, just plot scatter
        ax4.plot(rots, max_rate_rot.values, 'xk', mew=2, ms=10)
        ax4.plot(rots, max_rate_rot_5.values, 'xr', mew=2, ms=10)
        ax4.plot(rots, max_rate_rot_15.values, 'xb', mew=2, ms=10)
        ax5.plot(rots, max_rate_lat_rot.values, 'xk', mew=2, ms=10, label='10m')
        ax5.plot(rots, max_rate_lat_rot_5.values, 'xr', mew=2, ms=10, label='5m')
        ax5.plot(rots, max_rate_lat_rot_15.values, 'xb', mew=2, ms=10, label='15m')
    
    
    # Add approximate inverse relation to max rate plot, and approximations to trajectory plots
    f_crit = max_rate_lat_rot[2].values * rots[2]
    lat_predictions = f_crit/rots
    ax5.plot(rots, lat_predictions,'k')
    
    colors=['b','g','k','r','c','m','y']    
    for i in range(7):
        ax1.plot([lat_predictions[i],lat_predictions[i]], [-1.,1.], colors[i]+'--', alpha=0.5)
        ax2.plot([lat_predictions[i],lat_predictions[i]], [-1.,1.], colors[i]+'--', alpha=0.5)
        ax3.plot([lat_predictions[i],lat_predictions[i]], [-1.,1.], colors[i]+'--', alpha=0.5)
    
    # Add labels 
    ax4.set_ylabel('ITCZ migration rate')
    ax4.set_xlabel('$\Omega$/$\Omega_{E}$')
    ax4.set_yticks(np.arange(0.1,0.8,0.1))
    ax5.set_ylabel('Latitude')
    ax5.set_xlabel('$\Omega$/$\Omega_{E}$')
    ax5.set_yticks(np.arange(2,17,2))
    ax5.set_ylim([2,16])
    
    # Add grid
    ax4.grid(True,linestyle=':')
    ax5.grid(True,linestyle=':')
    
    legend = ax5.legend( loc='upper right', borderaxespad=0., fontsize=10, title='Mixed layer depth', ncol=2)
    legend.get_title().set_fontsize(10)
    
    # Add figure numbering
    ax1.text(-45, 1., 'a)')
    ax2.text(-45, 1., 'b)')
    ax3.text(-45, 1., 'c)')
    ax4.text(0.1, 0.7, 'd)')
    ax5.text(0.1, 16., 'e)')

    
    plt.subplots_adjust(left=0.075, right=0.97, top=0.95, bottom=0.1, hspace=0.4, wspace=1.)
    
    # Save as a pdf
    plt.savefig(plot_dir + 'pcent_rate_rot.pdf', format='pdf')
    plt.close()
    
    