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


def p_cent_rate(data, days=False):
    
    # Locate precipitation centroid
    precip_centroid(data)
    
    # If units in are days rather than pentads (e.g. for shorter years) convert to pentads
    if days:
        data['xofyear'] = data['xofyear']/5.
    
    # Get rate of movement of precip centroid
    dpcentdt = gr.ddt(data.p_cent) * 86400.
    
    dpcentdt2 = gr.ddt(dpcentdt, secperunit = 86400.) * 86400.
    
    return dpcentdt, dpcentdt2


def p_cent_rate_max(runs, days=None):
    # Get the maximum rate of change of the precipitation centroid latitude, and the latitude at which this occurs.
    
    if days==None:
        days = [False]*len(runs)
    i=0
    max_rate = []
    max_rate_lat = []
    
    for run in runs:
        # Open dataset
        data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
        dpcentdt, dpcentdt2 = p_cent_rate(data, days[i]) # Get rate of movement of precip centroid
    
        dpcentdt_ppos = dpcentdt.where(data.p_cent>=0.)   # Find precip centroid rate where precip centroid is in the northern hemisphere
        dpcentdt_max = dpcentdt_ppos.where(dpcentdt_ppos==dpcentdt_ppos.max(),drop=True)   # Find the maximum of the above
        pcent_dtmax = data.p_cent.sel(xofyear=dpcentdt_max.xofyear)    # Find the location of the preciptiation when the rate is maximum
        print dpcentdt_max.values, pcent_dtmax.values     # Print both
        
        max_rate.append(dpcentdt_max[0])
        max_rate_lat.append(pcent_dtmax[0])
    
    return max_rate, max_rate_lat
    
   
def p_cent_rate_plot(run, ax_in, ylim_in=0.55, days=False):
    # Evaluate the precipitation centroid rate and plot this on same axis as its location
    
    # Open dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    dpcentdt, dpcentdt2 = p_cent_rate(data, days) # Get rate of movement of precip centroid
    
    ax_twin = ax_in.twinx()  # Create a twin of the axes passed in

    data.p_cent.plot(ax=ax_twin, color='b')  # Plot the precip centroid location in blue on the twinned axes
    dpcentdt.plot(ax=ax_in, color='k')   # Plot the rate in black on the main axes
    
    # Set labels and limits
    ax_in.grid(True,linestyle=':')
    ax_in.set_xlabel('')
    ax_in.set_ylabel('')
    ax_twin.set_ylabel('')
    ax_twin.set_ylim([-27.5,27.5])
    ax_in.set_ylim([-1.*ylim_in, ylim_in])
    ax_in.set_title(run, fontsize=14)
    ax_twin.spines['right'].set_color('blue')
    plt.tight_layout()
    
    return ax_twin # Return twinned axes



def p_cent_grad_scatter(run, ax_in, ylim_in=0.55, days=False, period_fac=1.):
    
    # Open dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    dpcentdt, dpcentdt2 = p_cent_rate(data, days) # Get rate of movement of precip centroid
    
    ax_in.plot(data.p_cent, dpcentdt*period_fac, 'xk') # Plot scatter of precip centroid lat vs rate
    
    # Set labels and limits
    ax_in.grid(True,linestyle=':')
    ax_in.set_xlabel('')
    ax_in.set_ylabel('')
    ax_in.set_xlim([-25,25])
    ax_in.set_ylim([-1.*ylim_in, ylim_in])
    ax_in.set_title(run, fontsize=14)
    

def p_cent_acc_scatter(run, ax_in, ylim_in=0.5, days=False):
    
    # Open dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    dpcentdt, dpcentdt2 = p_cent_rate(data, days) # Get rate of movement of precip centroid
    
    ax_in.plot(data.p_cent, dpcentdt2, 'xk') # Plot scatter of precip centroid lat vs rate
    
    # Set labels and limits
    ax_in.grid(True,linestyle=':')
    ax_in.set_xlabel('')
    ax_in.set_ylabel('')
    ax_in.set_xlim([-25,25])
    ax_in.set_ylim([-1.*ylim_in, ylim_in])
    ax_in.set_title(run, fontsize=14)
    
    

if __name__ == "__main__":
    
    # Set plotting directory
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)

    rcParams['figure.figsize'] = 10, 5
    rcParams['font.size'] = 14
    
    
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
    axes = [ax1,ax2,ax3,ax4,ax5,ax6]
    runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
            'sn_2.000', 'sn_4.000', 'sn_8.000' ]
            
    # Set up list for twinned axes and do plots
    axes_twin=[]
    ax_twin = p_cent_rate_plot(runs[0], axes[0], days=True) 
    axes_twin.append(ax_twin)
    for i in range(1,6):
        ax_twin = p_cent_rate_plot(runs[i], axes[i])
        axes_twin.append(ax_twin)
    
    # Set x limits for axes
    axes[0].set_xlim([0,18])
    axes[1].set_xlim([0,36])
    axes[2].set_xlim([0,72])
    axes[3].set_xlim([0,144])
    axes[4].set_xlim([0,288])
    axes[5].set_xlim([0,576])
    
    # Remove clutter from axes labels
    for i in [0,1,3,4]:
        axes_twin[i].set_yticklabels('')
    for i in [1,2,4,5]:
        axes[i].set_yticklabels('')
    
    # Label axes where necessary
    for i in [0,3]:
        axes[i].set_ylabel('P. cent. lat. rate')
        
    for i in [2,5]:
        axes_twin[i].set_ylabel('Precip centroid lat.')
    
    for i in range(3,6):
        axes[i].set_xlabel('Pentad')
    
    for i in range(6):
        axes[i].tick_params(axis='both', which='major', labelsize=11)
        axes_twin[i].tick_params(axis='both', which='major', labelsize=11)
    
    plt.subplots_adjust(right=0.9, left=0.1, top=0.95, bottom=0.1, hspace=0.4, wspace=0.1)
    
    plt.savefig(plot_dir + 'precip_centroid_rate_seasons.pdf', format='pdf')
    plt.close()
    
    
    
    f2, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
    axes = [ax1,ax2,ax3,ax4,ax5,ax6]
    runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
            'sn_2.000', 'sn_4.000', 'sn_8.000' ]
    
    # Do plots
    p_cent_grad_scatter(runs[0], axes[0], days=True)
    for i in range(1,6):
        p_cent_grad_scatter(runs[i], axes[i])
    
    # Remove clutter from axes labels
    for i in [1,2,4,5]:
        axes[i].set_yticklabels('')
        
    # Label axes where necessary
    for i in range(0,6,3):
        axes[i].set_ylabel('P. cent. lat. rate')
    for i in range(3,6):
        axes[i].set_xlabel('Precip centroid lat.')
    
    for i in range(6):
        axes[i].tick_params(axis='both', which='major', labelsize=11)
        
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.1, hspace=0.4, wspace=0.2)
    
    plt.savefig(plot_dir + 'pcent_grad_scatter_seasons.pdf', format='pdf')
    plt.close()
    
    
    
    
    f2, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
    axes = [ax1,ax2,ax3,ax4,ax5,ax6]
    runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
            'sn_2.000', 'sn_4.000', 'sn_8.000' ]
    
    # Do plots
    period_fac = [0.25, 0.5, 1., 2., 4., 8.]
    p_cent_grad_scatter(runs[0], axes[0], days=True, period_fac=period_fac[0])
    for i in range(1,6):
        p_cent_grad_scatter(runs[i], axes[i], period_fac=period_fac[i])
    
    # Remove clutter from axes labels
    for i in [1,2,4,5]:
        axes[i].set_yticklabels('')
        
    # Label axes where necessary
    for i in range(0,6,3):
        axes[i].set_ylabel('P. cent. lat. rate')
    for i in range(3,6):
        axes[i].set_xlabel('Precip centroid lat.')
    
    for i in range(6):
        axes[i].tick_params(axis='both', which='major', labelsize=11)
        
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.1, hspace=0.4, wspace=0.2)
    
    plt.savefig(plot_dir + 'pcent_grad_scatter_seasons_scaled.pdf', format='pdf')
    plt.close()
    
    
    
    
    
    f2, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
    axes = [ax1,ax2,ax3,ax4,ax5,ax6]
    runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
            'sn_2.000', 'sn_4.000', 'sn_8.000' ]
    
    # Do plots
    p_cent_acc_scatter(runs[0], axes[0], days=True)
    for i in range(1,6):
        p_cent_acc_scatter(runs[i], axes[i])
    
    # Remove clutter from axes labels
    for i in [1,2,4,5]:
        axes[i].set_yticklabels('')
        
    # Label axes where necessary
    for i in range(0,6,3):
        axes[i].set_ylabel('P. cent. lat. acc.')
    for i in range(3,6):
        axes[i].set_xlabel('Precip centroid lat.')
    
    for i in range(6):
        axes[i].tick_params(axis='both', which='major', labelsize=11)
        
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.1, hspace=0.4, wspace=0.2)
    
    plt.savefig(plot_dir + 'pcent_acc_scatter_seasons.pdf', format='pdf')
    plt.close()
    
    
    
    
    
    
    

    rcParams['figure.figsize'] = 5, 7
    rcParams['font.size'] = 14
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    max_rate, max_rate_lat = p_cent_rate_max(runs)
    
    ax1.plot([0.25, 0.5, 1.0, 2.0, 4.0, 8.0], max_rate, 'xk')
    ax1.set_xlabel('')
    ax1.set_ylabel('Max rate')
    
    ax2.plot([0.25, 0.5, 1.0, 2.0, 4.0, 8.0], max_rate_lat, 'xk')
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')
    ax2.set_xlabel('Orbital period factor')
    ax2.set_ylabel('Lat of max rate')
    
    plt.subplots_adjust(right=0.95, left=0.15, top=0.95, bottom=0.1, hspace=0.1, wspace=0.2)
    
    plt.savefig(plot_dir + 'sn_scatter.pdf', format='pdf')
    plt.close()
    
    