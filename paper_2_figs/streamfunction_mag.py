"""
Make lat-time plots of the streamfunction magnitude for mld and rot runs 20/11/2018
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, model_constants as mc, make_sym
from climatology import precip_centroid
from hadley_cell import get_edge_psi
from pylab import rcParams
from hadley_cell import mass_streamfunction
from pcent_rate_max import p_cent_rate_max


rcParams['figure.figsize'] = 20, 5
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

def psi_plot(run, ax, title=''):    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    mse = (mc.cp_air * data.temp + mc.L * data.sphum + mc.grav * data.height).mean('lon').sel(pfull=850.)/1000.
    dmsedy = gr.ddy(mse, vector=False)*100000.
    dtsurfdy = gr.ddy(data.t_surf.mean('lon'), vector=False) * 100000.
    psi /= 1.e9
    psi = np.abs(psi).max('pfull')
    
    #mse.plot.contour(ax=ax, x='xofyear', y='lat', add_labels=False, levels=np.arange(200.,401.,10.), colors='w')
    dmsedy.plot.contourf(ax=ax, x='xofyear', y='lat', add_labels=False, levels=np.arange(-3.,3.1,0.25), add_colorbar=False)
    #dtsurfdy.plot.contourf(ax=ax, x='xofyear', y='lat', add_labels=False, levels=np.arange(-1.5,1.6,0.1), add_colorbar=False)
    psi.plot.contour(ax=ax, x='xofyear', y='lat', extend='both', add_labels=False, levels=np.arange(0,101.,100.), add_colorbar=False, alpha=0.4)
    ax.set_title(title, fontsize=17)
    ax.set_ylim(-60,60)
    ax.grid(True,linestyle=':')
    ax.set_yticks(np.arange(-60.,61.,30.))
    ax.set_xticks(np.arange(0.,72.,18.))


# Six subplots
fig, ((ax1, ax2, ax3, ax4, ax5, ax6), (ax7, ax8, ax9, ax10, ax11, ax12)) = plt.subplots(2, 6, sharex='col', sharey='row')
axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
plt.set_cmap('RdBu_r')
#First plot
psi_plot('rt_0.500', ax=ax1, title='rt0.5')
psi_plot('rt_0.750', ax=ax2, title='rt0.75')
psi_plot('rt_1.250', ax=ax3, title='rt1.25')
psi_plot('rt_1.500', ax=ax4, title='rt1.5')
psi_plot('rt_1.750', ax=ax5, title='rt1.75')
psi_plot('rt_2.000', ax=ax6, title='rt2')

psi_plot('mld_2.5', ax=ax7, title='mld2.5')
psi_plot('mld_5', ax=ax8, title='mld5')
f1 = psi_plot('sn_1.000', ax=ax9, title='Control')
psi_plot('mld_15', ax=ax10, title='mld15')
psi_plot('mld_20', ax=ax11, title='mld20')

ax1.set_ylabel('Latitude')
ax7.set_ylabel('Latitude')

for ax in [ax7,ax8,ax9,ax10,ax11]:
    ax.set_xlabel('Pentad')
    
#ax1.text(-15, 60, 'a)')
    
plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.12, hspace=0.25, wspace=0.12)
    #Colorbar
#cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    
plt.savefig(plot_dir + 'streamfunction_mag.pdf', format='pdf')
plt.close()

