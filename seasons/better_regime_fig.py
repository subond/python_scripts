"""
Plot overturning at 500 hPa at the Equator vs the interpolated precipitation centroid
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import numpy as np
import xarray as xr
#import sh
from data_handling import cell_area
#from pylab import rcParams
from physics import precip_centroid, model_constants as mc, mass_streamfunction, get_edge_psi
import scipy.interpolate as spint
from pylab import rcParams

month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 20
rcParams['text.usetex'] = True

# Get data for steady state experiments
intdown=True
vars_90 = get_edge_psi('ss_90.000', thresh=0., intdown=intdown)
vars_91 = get_edge_psi('ss_91.000', thresh=0., intdown=intdown)
vars_92 = get_edge_psi('ss_92.000', thresh=0., intdown=intdown)
vars_95 = get_edge_psi('ss_95.000', thresh=0., intdown=intdown)
vars_100 = get_edge_psi('ss_100.000', thresh=0., intdown=intdown)
vars_105 = get_edge_psi('ss_105.000', thresh=0., intdown=intdown)


def plot_regime(vars, varcolor='k', shem=True):
    if shem:
        psi_max = vars[1]*-1.
        edge_loc = vars[0]
    else:
        psi_max = vars[1]
        edge_loc = vars[0]*-1.

    i = int(vars[0].argmax('xofyear'))
    j = int(vars[0].argmin('xofyear'))
    plt.plot(edge_loc[j:i], psi_max[j:i], 'x', color=varcolor, ms=10, mew=2) #, alpha=0.7)
    plt.plot(edge_loc[i:], psi_max[i:], '+', color=varcolor, ms=10, mew=2, alpha=0.5)
    plt.plot(edge_loc[:j], psi_max[:j], '+', color=varcolor, ms=10, mew=2, alpha=0.5)


def plot3_regime(vars1, vars2, vars3, plotname, shem=True, logplot=True, add_ss=True):
    """Plot 3 sets of psi edge loc vs psi max on the same axis for comparison"""
    
    plot_regime(vars1, shem=shem)
    plot_regime(vars2, varcolor='b', shem=shem)
    plot_regime(vars3, varcolor='r', shem=shem)
    
    if add_ss:
        plt.plot(vars_90[0], -1.*vars_90[1], 'sk', ms=10)
        plt.plot(vars_91[0], -1.*vars_91[1], 'sk', ms=10)
        plt.plot(vars_92[0], -1.*vars_92[1], 'sk', ms=10)
        plt.plot(vars_95[0], -1.*vars_95[1], 'sk', ms=10)
        plt.plot(vars_100[0], -1.*vars_100[1], 'sk', ms=10)
        plt.plot(vars_105[0], -1.*vars_105[1], 'sk', ms=10)
    
    plt.xlabel('Cell edge')
    plt.ylabel('Max 500 hPa Mass Streamfunction')
    plt.grid(True,linestyle=':')
    if logplot:
        plt.xscale('log')
        plt.yscale('log')
        plt.minorticks_off()
        plt.xlim(0.625,40)
        plt.ylim(100,600)
        plt.xticks([1.25,2.5,5,10,20,40])
        plt.yticks([75,150,300,600])
        ax=plt.gca()
        ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(tk.ScalarFormatter())
        plt.tight_layout()
    else:
        plt.ylim([0,600])
        plt.xlim([0,40])
    plt.savefig('/scratch/rg419/plots/seasons/regimefig_' + plotname +'.pdf', format='pdf')
    plt.close()
            
            
vars_sn10 = get_edge_psi('sn_1.000')
vars_sn05 = get_edge_psi('sn_0.500')
vars_sn20 = get_edge_psi('sn_2.000')




plot3_regime(vars_sn10, vars_sn05, vars_sn20, 'seasons')


