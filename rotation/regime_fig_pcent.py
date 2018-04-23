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
from physics import precip_centroid, model_constants as mc, mass_streamfunction, psi_max_500
import scipy.interpolate as spint
from pylab import rcParams

month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 20
rcParams['text.usetex'] = True

# Get data for steady state experiments
intdown=True
ss_thresh = 0. 


def plot_regime(psi_max, pcent, varcolor='k', shem=True, guide=10):
    if shem:
        psi_max = psi_max*-1.
        edge_loc = pcent
    else:
        psi_max = psi_max
        edge_loc = pcent*-1.
    
    pcent = pcent.sel(xofyear=psi_max.xofyear)

    i = int(pcent[guide:].argmax('xofyear'))+guide
    j = int(pcent.argmin('xofyear'))
    #print j, i
    plt.plot(pcent[j:i], psi_max[j:i], 'x', color=varcolor, ms=10, mew=2) #, alpha=0.7)
    plt.plot(pcent[i:], psi_max[i:], '+', color=varcolor, ms=10, mew=2, alpha=0.5)
    plt.plot(pcent[:j], psi_max[:j], '+', color=varcolor, ms=10, mew=2, alpha=0.5)


def plot3_regime(vars_list, pcent_list, plotname, shem=True, logplot=True, add_ss=False, guide3=10):
    """Plot 3 sets of psi edge loc vs psi max on the same axis for comparison"""
    
    plot_regime(psimax_list[0], pcent_list[0], shem=shem)
    plot_regime(psimax_list[1], pcent_list[1], varcolor='b', shem=shem)
    plot_regime(psimax_list[2], pcent_list[2], varcolor='r', shem=shem, guide=guide3)
    
    if add_ss:
        plt.plot(vars_90[0], -1.*vars_90[1], 'sk', ms=10)
        plt.plot(vars_91[0], -1.*vars_91[1], 'sk', ms=10)
        plt.plot(vars_92[0], -1.*vars_92[1], 'sk', ms=10)
        plt.plot(vars_95[0], -1.*vars_95[1], 'sk', ms=10)
        plt.plot(vars_100[0], -1.*vars_100[1], 'sk', ms=10)
        plt.plot(vars_105[0], -1.*vars_105[1], 'sk', ms=10)
    
    plt.xlabel('Precipitation Centroid')
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
    plt.savefig('/scratch/rg419/plots/rotation/regimefig_pcent_' + plotname +'.pdf', format='pdf')
    plt.close()
            

psimax_list = []
pcent_list = []
#for run in ['sn_1.000', 'rt_0.500', 'rt_2.000']:
#for run in ['sn_1.000', 'sn_0.500', 'sn_2.000']:
for run in ['sn_1.000', 'ap_20', 'ap_2']:
#for run in ['ap_2', 'ap_20']:
    psimax_list.append(psi_max_500(run))
    pcent_list.append(precip_centroid(run))

#psimax_list.append(psi_max_500('full_qflux', lonin=[60.,150.]))
#pcent_list.append(precip_centroid('full_qflux', lonin=[60.,150.]))



plot3_regime(psimax_list, pcent_list, 'mld', add_ss=False)


