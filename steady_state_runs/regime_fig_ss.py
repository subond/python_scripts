# Load in data for steady state simulations and save average of years 11-15 in climatology folder. 

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from physics import get_edge_psi

intdown=True

vars_90 = get_edge_psi('ss_90.000', thresh=0., intdown=intdown)
vars_95 = get_edge_psi('ss_95.000', thresh=0., intdown=intdown)
vars_100 = get_edge_psi('ss_100.000', thresh=0., intdown=intdown)
vars_105 = get_edge_psi('ss_105.000', thresh=0., intdown=intdown)
vars_110 = get_edge_psi('ss_110.000', thresh=0., intdown=intdown)
vars_115 = get_edge_psi('ss_115.000', thresh=0., intdown=intdown)

plt.plot(vars_90[0], vars_90[1], 'xk', mew=2)
plt.plot(vars_95[0], vars_95[1], 'xk', mew=2)
plt.plot(vars_100[0], vars_100[1], 'xk', mew=2)
plt.plot(vars_105[0], vars_105[1], 'xk', mew=2)
plt.plot(vars_110[0], vars_110[1], 'xk', mew=2)
plt.plot(vars_115[0], vars_115[1], 'xk', mew=2)
plt.show()