''' 
16/08/2018 Evaluate Walker cell over a given latitude range
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import walker_cell
from data_handling_updates import model_constants as mc
from windspharm.xarray import VectorWind


def walker_cell_monthly(run, latin=None):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
        
    plot_dir = '/scratch/rg419/plots/overturning_monthly/' + run + '/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data.coords['month'] = (data.xofyear - 1) //6 + 1 
    data = data.groupby('month').mean(('xofyear'))
    
    # Create a VectorWind instance to handle the computation, and compute variables
    w = VectorWind(data.ucomp.sel(pfull=np.arange(50.,950.,50.)), data.vcomp.sel(pfull=np.arange(50.,950.,50.)))
    uchi, vchi, upsi, vpsi = w.helmholtz()
    
    psi_w = walker_cell(uchi, latin=latin, dp_in=-50.)
    psi_w /= 1.e9
    
    # Set figure parameters
    rcParams['figure.figsize'] = 12, 7
    rcParams['font.size'] = 14
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4)
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        psi_w.sel(month=i+1).plot.contour(ax=ax, x='lon', y='pfull', yincrease=False, levels=np.arange(0.,301,50.), colors='k', add_labels=False)
        psi_w.sel(month=i+1).plot.contour(ax=ax, x='lon', y='pfull', yincrease=False, levels=np.arange(-300.,0.,50.), colors='k', linestyles='dashed', add_labels=False)
        i=i+1
        ax.set_xticks(np.arange(0.,361.,120.))
        ax.set_yticks(np.arange(0.,1001.,250.))
        ax.grid(True,linestyle=':')
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    plt.savefig(plot_dir + 'walker_' + run + '.pdf', format='pdf')    
    plt.close()
    
walker_cell_monthly('q_shallow')
walker_cell_monthly('3q_shallow')
walker_cell_monthly('half_shallow')
walker_cell_monthly('half_shallow_5')
walker_cell_monthly('half_shallow_10')
