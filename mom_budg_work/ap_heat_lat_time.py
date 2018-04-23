# Produce plots of the zonal mean terms in the momentum budget as a function of time for the aquaplanet simulations, showing the transition and change in balance of terms

from physics import heatbudg_zon_pd_fn
from pentads import pentad_dic
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
plt.rc('text', usetex=True)

def plot_pd_heat_fn(inp_fol, years):
    
    data=heatbudg_zon_pd_fn(inp_fol, years,[0,72])
    pd_dic = pentad_dic(1)
    #data = rundata.mean(('lon'))
    
    def plot_heat_var(var,lev,levels):
        var_dic = {'heat_eddy': 'Eddy advective terms',
                    'heat_mean': 'Mean state advective terms'}
        plot_data = data.data_vars[var]
        plot_data = plot_data*86400.
        g=plot_data[:,lev,:].plot.contourf(x='pentad', y='lat', add_label = False, add_colorbar=False)
        plt.colorbar(g)
        tickspace = range(10,80,10)
        labels = [pd_dic[k] for k in tickspace]
        plt.ylim(0,90)
        plt.xlim(1,73)
        plt.xlabel('Pentad')
        plt.xticks(tickspace,labels,rotation=25)
        plt.ylabel('Latitude')
        plt.plot([float(data.pentad[36]),float(data.pentad[36])],[0.,90.],'k--')
        plt.title(var_dic[var])
        plt.savefig('/scratch/rg419/plots/mom_budg_work/'+inp_fol+'/'+var+'.png')
        plt.close()
        
    vars = ['heat_eddy','heat_mean']
    for var in vars:
        plot_heat_var(var,2,np.arange(-2.,2.1,0.2))            


plot_pd_heat_fn('aquaplanet_10m',range(21,41))

#plot_pd_mom_fn('aquaplanet_10m',range(21,41))