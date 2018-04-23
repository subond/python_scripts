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


data = time_means('flat_10m',[121,145],'pentad')
mombudg = mom_budg.mombudg_closure_fn(data)

