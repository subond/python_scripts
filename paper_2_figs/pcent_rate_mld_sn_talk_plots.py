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
    
    
def set_plot_features(ax, title='', legend_labels=[], fontsize=18, leg_title=None):
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
    rcParams['figure.figsize'] = 10, 7
    rcParams['font.size'] = 18

    # Set colors for lines
    colors=['C0','C1','C2','C3','C4','C5','C6']    
    
    runs_sn = ['sn_1.000', 'sn_2.000', 'sn_3.000', 'sn_4.000',
               'sn_0.500', 'sn_0.250']
    days = [False, False, False, False, False, True]
    period_fac = [1., 2., 3., 4., 0.5, 0.25]
    leg_labels = ['1.0', '2.0', '3.0', '4.0','0.5','0.25']
    
    
    subsolar_point = -23.439*np.cos(np.arange(0.,361.,1.)*np.pi/180.)
    subsolar_point = xr.DataArray(subsolar_point, coords=[np.arange(0.,361.,1.)], dims=['xofyear'])
    subsolar_rate = gr.ddt(subsolar_point, secperunit=86400.) * 86400.
    
    f, ax1 = plt.subplots()
    p_cent_grad_scatter('sn_1.000', color='k', ax=ax1, linewidth=4.)    
    ax1.plot(subsolar_point,subsolar_rate, color='m', linewidth=4.)
    set_plot_features(ax1, title='Varying Year Lengths', legend_labels=['Control ITCZ','Subsolar point'], leg_title='')
    ax1.set_ylabel('ITCZ migration rate')
    ax1.set_xlabel('ITCZ latitude')
    plt.savefig(plot_dir + 'pcent_rate_talk_solar.pdf', format='pdf')
    plt.close()
    
    
    '''Seasons'''
    for i in range(6):
        j=0
        f, ax1 = plt.subplots()
        name=runs_sn[i]
        for run in runs_sn[:i+1]:
            if run == 'sn_1.000':
                p_cent_grad_scatter(run, color='k', ax=ax1, linewidth=4.)    
            elif run == 'sn_0.250':
                p_cent_grad_scatter(run, days=True, color=colors[j-1], ax=ax1, linewidth=2.)
            else:
                p_cent_grad_scatter(run, color=colors[j-1], ax=ax1, linewidth=2.) 
            j=j+1
        set_plot_features(ax1, title='Varying Year Lengths', legend_labels=leg_labels[:i+1], leg_title='P/P$_{E}$')
        ax1.set_ylabel('ITCZ migration rate')
        ax1.set_xlabel('ITCZ latitude')
        plt.savefig(plot_dir + 'pcent_rate_talk_sn_' + name + '.pdf', format='pdf')
        plt.close()
            
    
    '''Seasons scaled'''
    for i in range(6):
        j=0
        f, ax1 = plt.subplots()
        name=runs_sn[i]
        for run in runs_sn[:i+1]:
            if run == 'sn_1.000':
                p_cent_grad_scatter(run, color='k', ax=ax1, linewidth=4.)    
            elif run == 'sn_0.250':
                p_cent_grad_scatter(run, days=True, period_fac=0.25, color=colors[j-1], ax=ax1, linewidth=2.)
            else:
                p_cent_grad_scatter(run, color=colors[j-1], period_fac=period_fac[j], ax=ax1, linewidth=2.) 
            j=j+1
        set_plot_features(ax1, title='Varying Year Lengths', legend_labels=leg_labels[:i+1], leg_title='P/P$_{E}$')
        ax1.set_ylabel('ITCZ migration rate')
        ax1.set_xlabel('ITCZ latitude')
        plt.savefig(plot_dir + 'pcent_rate_talk_sn_scaled_' + name + '.pdf', format='pdf')
        plt.close()

    
    '''Mixed layer depths'''
    
    runs_mld = ['sn_1.000', 'mld_15', 'mld_20', 'mld_5', 'mld_2.5']
    leg_labels = ['10.', '15.', '20.', '5.','2.5']
    
    for i in range(5):
        j=0
        f, ax1 = plt.subplots()
        name=runs_mld[i]
        for run in runs_mld[:i+1]:
            if run == 'sn_1.000':
                p_cent_grad_scatter(run, color='k', ax=ax1, linewidth=4.)    
            else:
                p_cent_grad_scatter(run, color=colors[j-1], ax=ax1, linewidth=2.) 
            j=j+1
        set_plot_features(ax1, title='Varying MLDs', legend_labels=leg_labels[:i+1], leg_title='Mixed layer depth, m')
        ax1.set_ylabel('ITCZ migration rate')
        ax1.set_xlabel('ITCZ latitude')
        plt.savefig(plot_dir + 'pcent_rate_talk_mld_' + name + '.pdf', format='pdf')
        plt.close()
        
        
    '''Rotations'''
    
    runs_rot = ['sn_1.000', 'rt_0.500', 'rt_0.750', 
            'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']
    leg_labels = ['1.0','0.5', '0.75', '1.25', '1.5', '1.75', '2.']
    
    for i in [6]:
        j=0
        f, ax1 = plt.subplots()
        name=runs_rot[i]
        for run in runs_rot[:i+1]:
            if run == 'sn_1.000':
                p_cent_grad_scatter(run, color='k', ax=ax1, linewidth=4.)    
            else:
                p_cent_grad_scatter(run, color=colors[j-1], ax=ax1, linewidth=2.) 
            j=j+1
        set_plot_features(ax1, title='', legend_labels=leg_labels[:i+1], leg_title='$\Omega$/$\Omega_{E}$')
        ax1.set_ylabel('ITCZ migration rate')
        ax1.set_xlabel('ITCZ latitude')
        plt.savefig(plot_dir + 'pcent_rate_talk_rot_' + name + '.pdf', format='pdf')
        plt.close()
    
    
    runs_5 = ['mld_5', 'rt_0.500_5', 'rt_0.750_5', 
            'rt_1.250_5', 'rt_1.500_5', 'rt_1.750_5', 'rt_2.000_5']
    
    for i in [6]:
        j=0
        f, ax1 = plt.subplots()
        name=runs_5[i]
        for run in runs_5[:i+1]:
            if run == 'mld_5':
                p_cent_grad_scatter(run, color='k', ax=ax1, linewidth=4.)    
            else:
                p_cent_grad_scatter(run, color=colors[j-1], ax=ax1, linewidth=2.) 
            j=j+1
        set_plot_features(ax1, title='', legend_labels=leg_labels[:i+1], leg_title='$\Omega$/$\Omega_{E}$')
        ax1.set_ylabel('ITCZ migration rate')
        ax1.set_xlabel('ITCZ latitude')
        plt.savefig(plot_dir + 'pcent_rate_talk_rot5_' + name + '.pdf', format='pdf')
        plt.close()
    
    
    runs_15 = ['mld_15', 'rt_0.500_15', 'rt_0.750_15', 
            'rt_1.250_15', 'rt_1.500_15', 'rt_1.750_15', 'rt_2.000_15']
    
    for i in [6]:
        j=0
        f, ax1 = plt.subplots()
        name=runs_15[i]
        for run in runs_15[:i+1]:
            if run == 'mld_15':
                p_cent_grad_scatter(run, color='k', ax=ax1, linewidth=4.)    
            else:
                p_cent_grad_scatter(run, color=colors[j-1], ax=ax1, linewidth=2.) 
            j=j+1
        set_plot_features(ax1, title='', legend_labels=['0.5', '0.75', '1.0', '1.25', '1.5', '1.75', '2.'], leg_title='$\Omega$/$\Omega_{E}$')
        ax1.set_ylabel('ITCZ migration rate')
        ax1.set_xlabel('ITCZ latitude')
        plt.savefig(plot_dir + 'pcent_rate_talk_rot15_' + name + '.pdf', format='pdf')
        plt.close()
    
    