# Plot EKE vs MAPE for different runs

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from eke import eke, trop_height, latmean
from mape import mape
from pylab import rcParams
import sh


def mape_eke_plot(run):
    
    rcParams['figure.figsize'] = 10, 6
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/mape_eke/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    eke_n, eke_s, eke_t = eke(run)
    mape_n, mape_s, mape_t = mape(run)
    
    plt.plot(mape_n/1.e6, eke_n/1.e6, 'xk')
    plt.title('Northern Hemisphere')
    plt.xlabel('MAPE, MJm$^{-2}$')
    plt.ylabel('EKE, MJm$^{-2}$')
    plt.xlim([0,10])
    plt.ylim([0,1.2])
    figname = 'mape_eke_n_' + run + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()

    plt.plot(mape_s/1.e6, eke_s/1.e6, 'xk')
    plt.title('Southern Hemisphere')
    plt.xlabel('MAPE, MJm$^{-2}$')
    plt.ylabel('EKE, MJm$^{-2}$')
    plt.xlim([0,10])
    plt.ylim([0,1.2])
    figname = 'mape_eke_s_' + run + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    plt.plot(mape_t/1.e6, eke_t/1.e6, 'xk')
    plt.xlabel('MAPE, MJm$^{-2}$')
    plt.ylabel('EKE, MJm$^{-2}$')
    plt.xlim([0,10])
    plt.ylim([0,1.2])
    figname = 'mape_eke_' + run + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    

#mape_eke_plot('ap_2')
#mape_eke_plot('ap_20')
mape_eke_plot('sn_1.000')
#mape_eke_plot('full_qflux')