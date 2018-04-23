# Produce plots of the zonal mean terms in the momentum budget as a function of time for the aquaplanet simulations, showing the transition and change in balance of terms

from data_handling import time_means
from physics import mom_budg
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

def load_data(inp_fol):
    data = time_means(inp_fol,[121,481],'pentad')
    mombudg = mom_budg.mombudg_fn(data)
    mombudg.to_netcdf(path='/scratch/rg419/plots/monsoon_analysis/'+inp_fol+'/mom_stat_data.nc', mode='w')

load_data('flat_10m')
load_data('topo_10m')
load_data('aquamountain_10m')
load_data('aquaplanet_10m')
load_data('aquaplanet_2m')

#mombudg.mom_trans[:,9,:,:].plot.contourf(x='lon', y='lat', col='xofyear', col_wrap=4, add_label = False)
#plt.ylim(-30,60)
#plt.xlim(60,180)

#mombudg.mom_stat[:,9,:,:].plot.contourf(x='lon', y='lat', col='xofyear', col_wrap=4, add_label = False)
#plt.ylim(-30,60)
#plt.xlim(60,180)

#mombudg.mom_mean[:,9,:,:].plot.contourf(x='lon', y='lat', col='xofyear', col_wrap=4, add_label = False)
#plt.ylim(-30,60)
#plt.xlim(60,180)

#plt.show()



