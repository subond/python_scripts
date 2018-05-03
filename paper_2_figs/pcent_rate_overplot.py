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
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')

    
    dpcentdt, dpcentdt2 = p_cent_rate(data, days) # Get rate of movement of precip centroid
        
    if ax == None:
        plt.plot(data.p_cent, dpcentdt*period_fac, color, linewidth=linewidth)
        #plt.plot([data.p_cent[-1],data.p_cent[0]], [dpcentdt[-1]*period_fac,dpcentdt[0]*period_fac], color, linewidth=linewidth)
    else:
        ax.plot(data.p_cent, dpcentdt*period_fac, color, linewidth=linewidth) # Plot scatter of precip centroid lat vs rate
        #ax.plot([data.p_cent[-1],data.p_cent[0]], [dpcentdt[-1]*period_fac,dpcentdt[0]*period_fac], color, linewidth=linewidth)


def set_plot_features(ax, title='', legend_labels=[], fontsize=14, leg_title=None):
    """
    Inputs:
        ax - axis to modify
        title - title for axis (default empty string)
        legend_labels - label for legend (default empty list)
    
    """
    # Shrink current axis by 10%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax.legend(legend_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=fontsize, title=leg_title)
    ax.set_xlim([-25,25])
    ax.set_ylim([-1., 1.]) 
    ax.set_title(title)
    ax.grid(True,linestyle=':')



if __name__ == "__main__":
    
    # Set plotting directory
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)

    # Set figure parameters
    rcParams['figure.figsize'] = 15, 8
    rcParams['font.size'] = 16

    # Start figure with 4 subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

    # Set colors for lines
    colors=['b','g','r','c','m','y']


    # First subplot - year length
    runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
            'sn_2.000', 'sn_4.000']
    # Set scaling factors for rate
    period_fac = [0.25, 0.5, 1., 2., 4.]
    # Do 0.25 first as this is in days
    p_cent_grad_scatter(runs[0], days=True, color=colors[0], ax=ax1, linewidth=1.)
    # Do other plots, plot sn_1.000 run as thicker line
    j=1
    for i in range(1,5):
        if runs[i] == 'sn_1.000':
            p_cent_grad_scatter(runs[i], color='k', ax=ax1, linewidth=2.)    
        else:
            p_cent_grad_scatter(runs[i], color=colors[j], ax=ax1, linewidth=1.) 
            j=j+1
    set_plot_features(ax1, title='Varying Year Lengths', legend_labels=['0.25', '0.5', '1.0', '2.0', '4.0'], leg_title='P/P$_{E}$')
    
    
    # Second subplot - year length, scaled
    runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
            'sn_2.000', 'sn_4.000']
    # Set scaling factors for rate
    period_fac = [0.25, 0.5, 1., 2., 4.]
    # Do 0.25 first as this is in days
    p_cent_grad_scatter(runs[0], days=True, period_fac=period_fac[0], color=colors[0], ax=ax2, linewidth=1.)
    # Do other plots, plot sn_1.000 run as thicker line
    j=1
    for i in range(1,5):
        if runs[i] == 'sn_1.000':
            p_cent_grad_scatter(runs[i], color='k', ax=ax2, linewidth=2.)    
        else:
            p_cent_grad_scatter(runs[i], period_fac=period_fac[i], color=colors[j], ax=ax2, linewidth=1.)    
            j=j+1
    set_plot_features(ax2, title='Varying Year Lengths (Rates scaled)', legend_labels=['0.25', '0.5', '1.0', '2.0', '4.0'], leg_title='P/P$_{E}$')
    
    
    # Third subplot - mixed layer depths
    runs = ['mld_2.5', 'mld_5', 'sn_1.000', 
            'mld_15', 'ap_20']
    j=0
    for i in range(0,5):
        if runs[i] == 'sn_1.000':
            p_cent_grad_scatter(runs[i], color='k', ax=ax3, linewidth=2.)    
        else:
            p_cent_grad_scatter(runs[i], color=colors[j], ax=ax3, linewidth=1.)
            j=j+1
    set_plot_features(ax3, title='Varying MLDs', legend_labels=['2.5', '5.', '10.', '15.', '20.'], leg_title='MLD (m)')
    
    # Fourth subplot - rotation rates
    #runs = ['rt_0.500_15','rt_0.750_15', 'mld_15', 
    #         'rt_1.250_15','rt_1.500_15','rt_1.750_15','rt_2.000_15']
    runs = ['rt_0.500', 'rt_0.750', 'sn_1.000', 
            'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']
    j=0
    for i in range(0,7):
        #if runs[i] == 'mld_15':
        if runs[i] == 'sn_1.000':
            p_cent_grad_scatter(runs[i], color='k', ax=ax4, linewidth=2.)    
        else:
            p_cent_grad_scatter(runs[i], color=colors[j], ax=ax4, linewidth=1.)
            j=j+1
    set_plot_features(ax4, title='Varying Rotation Rates', legend_labels=['0.5', '0.75', '1.0', '1.25', '1.5', '1.75', '2.0'], leg_title='$\Omega$/$\Omega_{E}$')
    
    
    for ax in [ax1, ax3]:
        ax.set_ylabel('P. cent. lat. rate')
    
    for ax in [ax3, ax4]:
        ax.set_xlabel('Precip centroid lat.')
    
    
    plt.subplots_adjust(left=0.07, right=0.9, top=0.95, bottom=0.1, hspace=0.25, wspace=0.37)
    
    # Save as a pdf
    plt.savefig(plot_dir + 'pcent_grad_scatter_overplots.pdf', format='pdf')
    plt.close()
    
    