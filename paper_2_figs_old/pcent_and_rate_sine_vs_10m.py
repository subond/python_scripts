"""
Plot precip centroid latitude and migration rate for 10m and sine sst experiments 02/02/2018

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_centroid
from data_handling_updates import gradients as gr
from pylab import rcParams
from pcent_rate_overplot import p_cent_rate



def p_cent_rate_plot(run, ax_in, ylim_in=1., days=False):
    # Evaluate the precipitation centroid rate and plot this on same axis as its location
    
    # Open dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    dpcentdt, dpcentdt2 = p_cent_rate(data, days) # Get rate of movement of precip centroid
    
    ax_twin = ax_in.twinx()  # Create a twin of the axes passed in

    data.p_cent.plot(ax=ax_twin, color='b', linewidth=2)  # Plot the precip centroid location in blue on the twinned axes
    dpcentdt.plot(ax=ax_in, color='k', linewidth=2)   # Plot the rate in black on the main axes
    
    # Set labels and limits
    ax_in.grid(True,linestyle=':')
    ax_in.set_xlabel('')
    ax_in.set_ylabel('')
    ax_twin.set_ylabel('')
    ax_twin.set_ylim([-20.,20.])
    ax_in.set_ylim([-1.*ylim_in, ylim_in])
    ax_twin.spines['right'].set_color('blue')
    plt.tight_layout()
    
    return ax_twin # Return twinned axes


# Set plotting directory
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

rcParams['figure.figsize'] = 5, 5
rcParams['font.size'] = 14
    
    
f, (ax1, ax2) = plt.subplots(2, sharex=True)
            
# Set up list for twinned axes and do plots
axes_twin=[]
#ax1_twin = p_cent_rate_plot('sn_1_sst', ax1) 
ax1_twin = p_cent_rate_plot('sn_1.000', ax1) 
ax2_twin = p_cent_rate_plot('sine_sst_10m', ax2) 

ax2.set_xlim([0,72])

#ax1_twin.set_yticklabels('')
#ax2.set_yticklabels('')
    
ax1.set_ylabel('P. cent. lat. rate')  
ax2.set_ylabel('P. cent. lat. rate')  
ax1_twin.set_ylabel('Precip centroid lat.', color='b')
ax2_twin.set_ylabel('Precip centroid lat.', color='b')

ax2.set_xticks([12,24,36,48,60,72])
ax1_twin.set_yticks([-20.,-10.,0.,10.,20.])
ax2_twin.set_yticks([-20.,-10.,0.,10.,20.])

ax1_twin.tick_params(axis='y', colors='blue')
ax2_twin.tick_params(axis='y', colors='blue')

ax2.set_xlabel('Pentad')

#axes[i].tick_params(axis='both', which='major', labelsize=11)
#axes_twin[i].tick_params(axis='both', which='major', labelsize=11)

plt.subplots_adjust(right=0.8, left=0.2, top=0.95, bottom=0.1, hspace=0.2)

plt.savefig(plot_dir + 'pcent_and_rate_sine_vs_10m.pdf', format='pdf')
plt.close()
    
