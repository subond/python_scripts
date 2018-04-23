"""
Plot overturning at 500 hPa at the Equator vs the interpolated precipitation centroid
"""

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as tk
import numpy as np
import xarray as xr
from physics import get_edge_psi_nh as get_edge_psi
from pylab import rcParams
import statsmodels.api as sm


month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 20
rcParams['text.usetex'] = True


def fit_power_law(edge_ll, psi_min, maxmin):
    
    lat = (edge_ll <= maxmin[1]) & (edge_ll > maxmin[0])
    edge = np.log( edge_ll[lat] )
    psi_min = np.log( psi_min[lat] )
    
    A = np.array([ edge, np.ones(edge.shape) ])
    #consts = np.linalg.lstsq( A.T, psi_min )[0] # obtaining the parameters
    #resids = np.linalg.lstsq( A.T, psi_min )[1] # obtaining the parameters
    #print psi_min
    #print A


    model = sm.OLS(psi_min.values, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    
    #print result.summary()
        
    print '=== Coeffs ==='
    print np.exp(consts[1]), consts[0]
    print '=== Std Errs ==='
    print 2*std_err[1]*np.exp(consts[1]), 2*std_err[0]
    print '=== No obs ==='
    print result.nobs
    
    line = np.exp(consts[1]) * np.arange(maxmin[0],maxmin[1]+1)**(consts[0])
    
    return line, consts
    
    
    

def plot_regime(vars, varcolor='k', shem=True, guide=10, symbol='x', do_linefit_l=False, do_linefit_h=False):
    # shem specfies whether northern or southern hemisphere is being plotted
    # guide helps in locaton of outward and return branch, selected for shem=True
    if shem:
        psi_max = vars[1]#*-1.
        edge_loc = vars[0]
    else:
        psi_max = vars[1]
        edge_loc = vars[0]*-1.

    i1 = int(vars[0][guide:].argmax('xofyear'))+guide
    i2 = int(vars[1][guide:].argmax('xofyear'))+guide
    i = min([i1,i2])
    #i=i+1
    j = int(vars[0][0:50].argmin('xofyear'))
    
    print j, i
    plt.plot(edge_loc[j:i], psi_max[j:i], symbol, color=varcolor, ms=10, mew=2) #, alpha=0.7)
    if do_linefit_l:
        line, consts = fit_power_law(edge_loc[j:i], psi_max[j:i], maxmin=[0.5,7.])
        line = np.exp(consts[1]) * np.arange(0.5,7.5,0.5)**(consts[0])
        plt.plot(np.arange(0.5,7.5,0.5), line, varcolor+':', linewidth=2)
    if do_linefit_h:
        line, consts = fit_power_law(edge_loc[j:i], psi_max[j:i], maxmin=[7.,30.])
        plt.plot(np.arange(7.,31.), line, varcolor+':', linewidth=2)
    #plt.plot(edge_loc[i:], psi_max[i:], '+', color=varcolor, ms=10, mew=2, alpha=0.5)
    #plt.plot(edge_loc[:j], psi_max[:j], '+', color=varcolor, ms=10, mew=2, alpha=0.5)


def plot3_regime(vars1, vars2, vars3, era, vars4, plotname, shem=True, logplot=True, add_ss=True, guide3=10, key_labels=None):
    """Plot 3 sets of psi edge loc vs psi max on the same axis for comparison"""
    
    #plot_regime(vars1, shem=shem, do_linefit_h=True)
    #plot_regime(vars2, varcolor='b', shem=shem, do_linefit_l=True)
    #plot_regime(vars3, varcolor='r', shem=shem, guide=guide3, symbol='+', do_linefit_h=True, do_linefit_l=True)
    plot_regime(vars1, shem=shem, symbol='x')#, do_linefit_h=True)
    plot_regime(vars2, shem=shem, symbol='+')#, do_linefit_l=True)
    plot_regime(vars3, shem=shem, symbol='o', do_linefit_h=True, do_linefit_l=True)
    plot_regime(vars4, varcolor='r', shem=shem, guide=guide3, do_linefit_h=True, do_linefit_l=True)
    #plot_regime(era, varcolor='r', shem=shem, guide=guide3, symbol='s')
    
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
        plt.xlim(0.5,35)
        plt.ylim(150,4000)
        plt.xticks([1,5,10,20,30])
        plt.yticks([150, 200,300,400,500])
        ax=plt.gca()
        ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(tk.ScalarFormatter())
        plt.tight_layout()
    else:
        plt.ylim([0,600])
        plt.xlim([0,40])
    if not key_labels==None:
        plt.legend(handles=key_labels, loc='upper left', prop={'size': 18}, numpoints=1)
    
    plt.savefig('/scratch/rg419/plots/paper_1_figs/revisions/regimefig_' + plotname +'.pdf', format='pdf')
    plt.close()
            
thresh = 120.   
lev = 500. 
vars_ap2 = get_edge_psi('ap_2', thresh=thresh, lev=lev, sanity_check=False)
vars_ap10 = get_edge_psi('sn_1.000', thresh=thresh, lev=lev)
vars_ap20 = get_edge_psi('ap_20', thresh=thresh, lev=lev, sanity_check=False)
#vars_full_all = get_edge_psi('full_qflux', thresh=thresh, lev=lev)#, lonin=[60.,150.])
vars_full = get_edge_psi('full_qflux', thresh=thresh, lonin=[60.,150.], lev=lev)#, do_ageostrophic=True, sanity_check=True)
#vars_full = get_edge_psi('full_qflux', thresh=thresh, lonin=[345.,45.], lev=lev)
vars_era = get_edge_psi('era', lonin=[60.,150.], thresh=thresh, lev=lev, sanity_check=False)

    
key_labels = [mlines.Line2D(range(1), range(1), marker='o', color='white',
                      markersize=10, mew=2, markeredgecolor='black',
                      markerfacecolor='black', label='ap10')]
key_labels.append(mlines.Line2D(range(1), range(1), marker='x', color='white',
                      markersize=10, mew=2, markeredgecolor='black',
                      markerfacecolor='white', label='ap2'))
key_labels.append(mlines.Line2D(range(1), range(1), marker='+', color='white',
                      markersize=10, mew=2, markeredgecolor='black',
                      markerfacecolor='white', label='ap20'))
key_labels.append(mlines.Line2D(range(1), range(1), marker='x', color='white',
                      markersize=10, mew=2, markeredgecolor='red',
                      markerfacecolor='white', label='full'))
                          
plot3_regime(vars_ap2, vars_ap20, vars_ap10, vars_era, vars_full, 'ap_to_full', guide3=30, add_ss=False, key_labels=key_labels)





