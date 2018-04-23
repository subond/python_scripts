# Produce plots of the zonal mean terms in the momentum budget as a function of time for the aquaplanet simulations, showing the transition and change in balance of terms

from physics import mombudg_pd_fn
from plotting import load_bs_vars
from data_handling import pentad_dic, month_dic
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


plt.rc('text', usetex=True)
font = {'size'   : 18}
plt.rc('font', **font)
mathtxt={'fontset':'custom', 'default':'regular'}
plt.rc('mathtext',**mathtxt)


def plot_pd_mom_fn(inp_fol, years):
    data_p = load_bs_vars(inp_fol, years)
    data=mombudg_pd_fn(inp_fol, years)
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    def plot_mom_var(var,lev,levels):
        var_dic = {'fv_av': 'fv',
                    'mom_eddy': 'Eddy advective terms',
                    'mom_mean': 'Mean state advective terms',
                    'mom_sum': 'Residual'}
        plot_data = data.data_vars[var]
        plot_data = plot_data*10000.
        g=plot_data[:,lev,:].plot.contourf(x='pentad', y='lat',levels=levels, add_label = False, add_colorbar=False, extend='neither')
        cb1=plt.colorbar(g)
        cb1.set_label('$\displaystyle10^{-4}m/s^2$')
        cs=data_p.totp.plot.contour(x='pentad', y='lat',levels=np.arange(6.,31.,6.), add_label = False, add_colorbar=False,colors='k')
        plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
        plt.ylim(0,90)
        plt.xlim(1,73)
        plt.xlabel('')
        plt.xticks(tickspace,labels,rotation=25)
        plt.ylabel('Latitude')
        #plt.plot([float(data.pentad[36]),float(data.pentad[36])],[0.,90.],'k--')
        plt.tight_layout()
        plt.title('')
        #plt.title(var_dic[var] + ', $\displaystyle10^{-4}m/s^2$')
        plt.savefig('/scratch/rg419/plots/mom_budg_work/'+inp_fol+'/'+var+'_lattime.png')
        plt.close()
        
    vars = ['fv_av','mom_eddy','mom_mean','mom_sum']
    for var in vars:
        plot_mom_var(var,9,np.arange(-2.,2.1,0.2))            
        #plot_mom_var(var,9,np.arange(-0.0002,0.00021,0.00002))            

plot_pd_mom_fn('aquaplanet_2m',range(11,41))

#plot_pd_mom_fn('aquaplanet_10m',range(21,41))