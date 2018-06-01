"""
Produce plots showing precip centroid lat vs rate for experiments with varying obliquity

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
    #ax.set_ylim([-1., 1.]) 
    #ax.set_title(title, fontsize=16)
    ax.grid(True,linestyle=':')



if __name__ == "__main__":
    
    # Set plotting directory
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)

    # Set figure parameters
    rcParams['figure.figsize'] = 5.5, 10
    rcParams['font.size'] = 14

    # Start figure with 4 subplots
    fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1)

    # Set colors for lines
    colors=['b','g','r','c','m','y']
    
    runs = ['ob_10.000', 'ob_20.000', 'sn_1.000', 'ob_30.000', 'ob_40.000']
    obs = np.array([10.,20., 23.439, 30., 40.])
    
    i=0
    max_rate_ob, max_rate_lat_ob, max_lat_ob = p_cent_rate_max(runs)
    print('ob peak lat:', np.mean(max_rate_lat_ob[1:].values), '+-', np.std(max_rate_lat_ob[1:].values))
    dpdt_eq_ob = rate_at_eq(runs)
    for run in runs:
        if run == 'sn_1.000':
            p_cent_grad_scatter(run, color='k', ax=ax1, linewidth=2.)    
        else:
            p_cent_grad_scatter(run, color=colors[i], ax=ax1, linewidth=1.)
            i=i+1
    set_plot_features(ax1, title='', legend_labels=['10.','20.','23.439','30.', '40.'], leg_title='Obliquity')
    ax1.set_ylabel('ITCZ migration rate')
    ax1.set_xlabel('Precip centroid lat.')
    
    ax2.plot(obs, max_rate_ob.values, 'xk', mew=2, ms=10)
    #ax2.plot(obs, dpdt_eq_ob, '+k', mew=2, ms=10)
    ax2.set_ylabel('ITCZ migration rate')
    ax2.set_xlabel('Obliquity')
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')
    
    ax3.plot(obs, max_rate_lat_ob.values, 'xk', mew=2, ms=10)
    ax3.set_ylabel('Latitude')
    ax3.set_xlabel('Obliquity')
    #ax3.set_xscale('log')
    #ax3.set_yscale('log')
    
    ax2.grid(True,linestyle=':')
    ax3.grid(True,linestyle=':')
    ax3.set_yticks(np.arange(4,12,1))
    
    plt.subplots_adjust(left=0.15, right=0.9, top=0.95, bottom=0.07, hspace=0.3)
    
    # Save as a pdf
    plt.savefig(plot_dir + 'pcent_rate_ob.pdf', format='pdf')
    plt.close()
    
    