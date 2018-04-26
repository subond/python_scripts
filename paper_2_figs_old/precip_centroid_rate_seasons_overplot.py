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
from pcent_rate_overplot import p_cent_grad_scatter, set_plot_features

# Set plotting directory
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

rcParams['figure.figsize'] = 10, 5
rcParams['font.size'] = 20


   
# Set colors for lines
colors=['r','m','k','c','b']

fig = plt.figure()
ax1 = plt.subplot(111)

runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
        'sn_2.000', 'sn_4.000']
# Set scaling factors for rate
period_fac = [0.25, 0.5, 1., 2., 4.]
# Do 0.25 first as this is in days
p_cent_grad_scatter(runs[0], days=True, period_fac=period_fac[0], color=colors[0], ax=ax1, linewidth=2.)
# Do other plots, plot sn_1.000 run as thicker line
for i in range(1,5):
    if runs[i] == 'sn_1.000':
        p_cent_grad_scatter(runs[i], period_fac=period_fac[i], color=colors[i], ax=ax1, linewidth=4.)    
    else:
        p_cent_grad_scatter(runs[i], period_fac=period_fac[i], color=colors[i], ax=ax1, linewidth=2.)    
set_plot_features(ax1, title='', legend_labels=['0.25', '0.5', '1.0', '2.0', '4.0', '8.0'], fontsize=16)

plt.xlabel('Precip centroid lat.')

plt.subplots_adjust(left=0.15, right=0.8, top=0.95, bottom=0.15)

plt.savefig(plot_dir + 'seasons_overplot_scaled.pdf', format='pdf')
plt.close()
    

fig = plt.figure()
ax1 = plt.subplot(111)

runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
        'sn_2.000', 'sn_4.000']
# Set scaling factors for rate
period_fac = [0.25, 0.5, 1., 2., 4.]
# Do 0.25 first as this is in days
p_cent_grad_scatter(runs[0], days=True, color=colors[0], ax=ax1, linewidth=2.)
# Do other plots, plot sn_1.000 run as thicker line
for i in range(1,5):
    if runs[i] == 'sn_1.000':
        p_cent_grad_scatter(runs[i], color=colors[i], ax=ax1, linewidth=4.)    
    else:
        p_cent_grad_scatter(runs[i], color=colors[i], ax=ax1, linewidth=2.)    
set_plot_features(ax1, title='', legend_labels=['0.25', '0.5', '1.0', '2.0', '4.0', '8.0'], fontsize=16)

plt.xlabel('Precip centroid lat.')

plt.subplots_adjust(left=0.15, right=0.8, top=0.95, bottom=0.15)

plt.savefig(plot_dir + 'seasons_overplot.pdf', format='pdf')
plt.close()


rcParams['figure.figsize'] = 10, 10

fig, (ax1, ax2) = plt.subplots(2, sharex=True)

# Do 0.25 first as this is in days
p_cent_grad_scatter(runs[0], days=True, color=colors[0], ax=ax1, linewidth=2.)
# Do other plots, plot sn_1.000 run as thicker line
for i in range(1,5):
    if runs[i] == 'sn_1.000':
        p_cent_grad_scatter(runs[i], color=colors[i], ax=ax1, linewidth=4.)    
    else:
        p_cent_grad_scatter(runs[i], color=colors[i], ax=ax1, linewidth=2.)    
set_plot_features(ax1, title='', legend_labels=['0.25', '0.5', '1.0', '2.0', '4.0', '8.0'], fontsize=16)


# Do 0.25 first as this is in days
p_cent_grad_scatter(runs[0], days=True, period_fac=period_fac[0], color=colors[0], ax=ax2, linewidth=2.)
# Do other plots, plot sn_1.000 run as thicker line
for i in range(1,5):
    if runs[i] == 'sn_1.000':
        p_cent_grad_scatter(runs[i], period_fac=period_fac[i], color=colors[i], ax=ax2, linewidth=4.)    
    else:
        p_cent_grad_scatter(runs[i], period_fac=period_fac[i], color=colors[i], ax=ax2, linewidth=2.)    
set_plot_features(ax2, title='', legend_labels=['0.25', '0.5', '1.0', '2.0', '4.0', '8.0'], fontsize=16)


ax2.set_xlabel('Precip centroid lat.')

plt.subplots_adjust(left=0.15, right=0.8, top=0.95, bottom=0.1)

plt.savefig(plot_dir + 'seasons_overplot_both.pdf', format='pdf')
plt.close()

