"""
6/11/2017
Plot max overturning at 500 hPa for cross-equatorial cell vs the latitude of max near surface mse, or vs precip centroid, or vs lat at which cell drops below some threshold
Version for steady state runs
"""

import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import numpy as np
import xarray as xr
from climatology import peak_mse, precip_centroid
from pylab import rcParams
from hadley_cell import get_edge_psi
import statsmodels.api as sm


lev = 500. 
    
rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 20
rcParams['text.usetex'] = True


def set_vars(data, plottype=None, lonin=[-1.,361.], thresh=0.):
    # Load up variables to plot, depending on plot type
    edge_loc, psi_max, psi_max_loc = get_edge_psi(data, thresh=thresh, lev=lev, lonin=lonin, nh=True, sanity_check=True)
        
    if plottype == 'mse':
        data_mse = peak_mse(data, lonin=lonin)
        vars_out = (data_mse.mse_max_loc, psi_max)
        
    elif plottype == 'pcent':
        data = precip_centroid(data,lonin=lonin)
        vars_out = (data.p_cent, psi_max)

    elif plottype == None:
        vars_out = (edge_loc, psi_max)
    
    return vars_out
        

def load_vars(runs, plottype=None, lonin=[-1.,361.], thresh=0.):
    # load variables for multiple runs
    vars_out = []
    for run in runs:
        data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
        run_vars = set_vars(data, plottype=plottype, lonin=lonin, thresh=thresh)
        vars_out.append(run_vars)
    return vars_out


def plot_ss_regime(var_list, plotname, logplot=True, vc=None, s = None, key_labels=None):
    
    n = len(var_list)
    
    if s==None:
        s = ['x']*n
    if vc == None:
        vc = ['k']*n
    
    for i in range(n):
        plt.plot(var_list[i][0], var_list[i][1], color=vc[i], marker = s[i], markersize=10, mew=2.)
    
    plt.xlabel('Cell edge')
    plt.ylabel('Max 500 hPa Mass Streamfunction')
    plt.grid(True,linestyle=':')
    if logplot:
        plt.xscale('log')
        plt.yscale('log')
        plt.minorticks_off()
        #plt.xlim(0.05,40)
        #plt.ylim(100,500)
        plt.xlim(0.05,50)
        plt.ylim(50,900)
        plt.xticks([0.2,1,5,10,20,40])
        plt.yticks([100, 200,400,800])
        #plt.xlim(0.5,int(latmax)+5)
        #plt.ylim(int(psimin),int(psimax)+50)
        #xticklist = [i for i in [1,2,5,10,20,30,40,50] if i < latmax]
        #plt.xticks(xticklist)
        #plt.yticks(range(int(round(psimin/50)*50), int(round(psimax/50)*50), 100))
        ax=plt.gca()
        ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(tk.ScalarFormatter())
        plt.tight_layout()
    #else:
        #plt.ylim([0,int(latmax)+5])
        #plt.xlim([0,int(psimax)+50])
    if not key_labels==None:
        #handles, labels = ax.get_legend_handles_labels()
        #print handles
        plt.legend(handles=key_labels, loc='upper left', prop={'size': 12}, numpoints=1)
    plt.savefig('/scratch/rg419/plots/regime_fig/' + plotname +'.pdf', format='pdf')
    plt.close()






if __name__ == "__main__":
    
    #lonin = [-1.,361.]
    #runs = ['ss_90.000', 'ss_91.000', 'ss_92.000', 'ss_93.000', 'ss_94.000', 'ss_95.000', 'ss_100.000', 'ss_105.000', 'ss_110.000', 'ss_115.000']
    
    #vars_out = load_vars(runs)
    #plot_ss_regime(vars_out, 'ss_symmetric')
    
    #runs = ['qflux_0_100', 'qflux_s_100', 'qflux_10_100', 'qflux_15_100', 'qflux_20_100', 'qflux_25_100',
    #        'qflux_0_200', 'qflux_5_200', 'qflux_10_200', 'qflux_15_200', 'qflux_20_200', 'qflux_25_200',
    #        'qflux_0_300', 'qflux_5_300', 'qflux_10_300', 'qflux_15_300', 'qflux_20_300', 'qflux_25_300']
    #colors = ['b']*6 + ['y']*6 + ['r']*6
    #vars_out = load_vars(runs)
    #plot_ss_regime(vars_out, 'ss_qfluxes', vc=colors)
    
    
    runs = ['ss_90.000', 'ss_91.000', 'ss_92.000', 'ss_93.000', 'ss_94.000', 'ss_95.000', 'ss_100.000', 'ss_105.000', 'ss_110.000', 'ss_115.000',
            'qflux_0_100', 'qflux_5_100', 'qflux_10_100', 'qflux_15_100', 'qflux_20_100', 'qflux_25_100',
            'qflux_0_200', 'qflux_5_200', 'qflux_10_200', 'qflux_15_200', 'qflux_20_200', 'qflux_25_200',
            'qflux_0_300', 'qflux_5_300', 'qflux_10_300', 'qflux_15_300', 'qflux_20_300', 'qflux_25_300']
    
    #colors = ['k']*10 + ['b']*6 + ['y']*6 + ['r']*6
    colors = ['k']*10 + ['r','m','y','g','c','b']*3
    marker_type = ['x']*10 + ['+']*6 + ['o']*6 + ['s']*6
    vars_out = load_vars(runs, lonin=[275.,50.])
    
    key_labels = [mlines.Line2D(list(range(1)), list(range(1)), marker='x', color='white',
                          markersize=10, mew=2, markeredgecolor='black',
                          markerfacecolor='white', label='Symmetric runs')]
    key_labels.append(mpatches.Patch(color='m', label='5$^{\circ}$'))
    key_labels.append(mpatches.Patch(color='y', label='10$^{\circ}$'))
    key_labels.append(mpatches.Patch(color='g', label='15$^{\circ}$'))
    key_labels.append(mpatches.Patch(color='c', label='20$^{\circ}$'))
    key_labels.append(mpatches.Patch(color='b', label='25$^{\circ}$'))
    key_labels.append(mlines.Line2D(list(range(1)), list(range(1)), marker='+', color='white',
                          markersize=10, mew=2, markeredgecolor='black',
                          markerfacecolor='white', label='100 W/m$^2$'))
    key_labels.append(mlines.Line2D(list(range(1)), list(range(1)), marker='o', color='white',
                          markersize=10, mew=2, markeredgecolor='black',
                          markerfacecolor='white', label='200 W/m$^2$'))
    key_labels.append(mlines.Line2D(list(range(1)), list(range(1)), marker='s', color='white',
                          markersize=10, mew=2, markeredgecolor='black',
                          markerfacecolor='white', label='300 W/m$^2$'))
    
    
    plot_ss_regime(vars_out, 'ss_all_275_50', vc=colors, s=marker_type, key_labels=key_labels)




