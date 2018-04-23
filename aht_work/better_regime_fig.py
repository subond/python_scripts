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
run_list = ['ss_90.000', 'ss_91.000', 'ss_92.000', 'ss_93.000', 'ss_94.000', 'ss_95.000', 'ss_100.000', 'ss_105.000']
intdown=True
ss_thresh = 0. 

ss_vars=[]
for run in run_list:
    ss_vars.append(get_edge_psi(run, thresh=ss_thresh, intdown=intdown))
    

def plot_regime(vars, varcolor='k', shem=True, guide=10):
    # shem specfies whether northern or southern hemisphere is being plotted
    # guide helps in locaton of outward and return branch, selected for shem=True
    if shem:
        psi_max = vars[1]*-1.
        edge_loc = vars[0]
    else:
        psi_max = vars[1]
        edge_loc = vars[0]*-1.

    i = int(vars[0][guide:].argmax('xofyear'))+guide
    j = int(vars[0].argmin('xofyear'))
    #print j, i
    plt.plot(edge_loc[j:i], psi_max[j:i], 'x', color=varcolor, ms=10, mew=2) #, alpha=0.7)
    plt.plot(edge_loc[i:], psi_max[i:], '+', color=varcolor, ms=10, mew=2, alpha=0.5)
    plt.plot(edge_loc[:j], psi_max[:j], '+', color=varcolor, ms=10, mew=2, alpha=0.5)


def plot3_regime(vars1, vars2, vars3, plotname, shem=True, logplot=True, add_ss=True, guide3=10):
    """Plot 3 sets of psi edge loc vs psi max on the same axis for comparison"""
    
    plot_regime(vars1, shem=shem)
    plot_regime(vars2, varcolor='b', shem=shem)
    plot_regime(vars3, varcolor='r', shem=shem, guide=guide3)
    
    if add_ss:
        for i in range(len(ss_vars)):
            plt.plot(ss_vars[i][0], -1.*ss_vars[i][1], 'sk', ms=10)

    
    plt.xlabel('Cell edge')
    plt.ylabel('Max 500 hPa Mass Streamfunction')
    plt.grid(True,linestyle=':')
    if logplot:
        plt.xscale('log')
        plt.yscale('log')
        plt.minorticks_off()
        #plt.xlim(0.625,40)
        #plt.ylim(100,600)
        #plt.xticks([1.25,2.5,5,10,20,40])
        #plt.yticks([75,150,300,600])
        plt.xlim(1,30)
        plt.ylim(175,500)
        plt.xticks([1,5,10,20,30])
        plt.yticks([200,300,400,500])
        ax=plt.gca()
        ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(tk.ScalarFormatter())
        plt.tight_layout()
    else:
        plt.ylim([0,600])
        plt.xlim([0,40])
    plt.savefig('/scratch/rg419/plots/aht_work/regimefig_' + plotname +'.pdf', format='pdf')
    plt.close()
            
            
vars_ap2 = get_edge_psi('ap_2', thresh=120.)
vars_ap20 = get_edge_psi('ap_20', thresh=120.)
vars_full = get_edge_psi('full_qflux', lonin=[60.,150.], thresh=120.)

plot3_regime(vars_ap2, vars_ap20, vars_full, 'ap_to_full', guide3=30, add_ss=False)
#plot3_regime(vars_ap2, vars_ap20, vars_full, 'ap_to_full_lin', logplot=False)
#plot3_regime(vars_ap2, vars_ap20, vars_full, 'ap_to_full_nhem', shem=False)

#vars_ap2 = get_edge_psi('ap_2')
#vars_ap20 = get_edge_psi('ap_20')
#vars_sn10 = get_edge_psi('sn_1.000')
#vars_sn05 = get_edge_psi('sn_0.500')
#vars_sn20 = get_edge_psi('sn_2.000')
#vars_rt05 = get_edge_psi('rt_0.500')
#vars_rt20 = get_edge_psi('rt_2.000')

#plot3_regime(vars_ap2, vars_ap20, vars_sn10, 'mld')
#plot3_regime(vars_ap2, vars_ap20, vars_sn10, 'mld_lin', logplot=False)
#plot3_regime(vars_ap2, vars_ap20, vars_sn10, 'mld_nhem', shem=False)

#plot3_regime(vars_sn10, vars_sn05, vars_sn20, 'seasons')
#plot3_regime(vars_sn10, vars_sn05, vars_sn20, 'seasons_lin', logplot=False)
#plot3_regime(vars_sn10, vars_sn05, vars_sn20, 'seasons_nhem', shem=False)

#plot3_regime(vars_sn10, vars_rt05, vars_rt20, 'rotation', add_ss=False)
#plot3_regime(vars_sn10, vars_rt05, vars_rt20, 'rotation_lin', logplot=False)
#plot3_regime(vars_sn10, vars_rt05, vars_rt20, 'rotation_nhem', shem=False)

