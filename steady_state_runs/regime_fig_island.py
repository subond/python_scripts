# Load in data for steady state simulations and save average of years 11-15 in climatology folder. 

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from physics import get_edge_psi

intdown=True

lonin = [-1.,361.]

vars_0_100 = get_edge_psi('qflux_0_100', lonin=lonin, thresh=0., intdown=intdown)
vars_5_100 = get_edge_psi('qflux_5_100', lonin=lonin, thresh=0., intdown=intdown)
vars_10_100  = get_edge_psi('qflux_10_100', lonin=lonin, thresh=0., intdown=intdown)
vars_15_100  = get_edge_psi('qflux_15_100', lonin=lonin, thresh=0., intdown=intdown)
vars_20_100  = get_edge_psi('qflux_20_100', lonin=lonin, thresh=0., intdown=intdown)
vars_25_100  = get_edge_psi('qflux_25_100', lonin=lonin, thresh=0., intdown=intdown)

vars_0_200  = get_edge_psi('qflux_0_200', lonin=lonin, thresh=0., intdown=intdown) #
vars_5_200 = get_edge_psi('qflux_5_200', lonin=lonin, thresh=0., intdown=intdown) 
vars_10_200 = get_edge_psi('qflux_10_200', lonin=lonin, thresh=0., intdown=intdown)
vars_15_200 = get_edge_psi('qflux_15_200', lonin=lonin, thresh=0., intdown=intdown)
vars_20_200 = get_edge_psi('qflux_20_200', lonin=lonin, thresh=0., intdown=intdown)#
vars_25_200 = get_edge_psi('qflux_25_200', lonin=lonin, thresh=0., intdown=intdown)#

vars_0_300 = get_edge_psi('qflux_0_300', lonin=lonin, thresh=0., intdown=intdown)#
vars_5_300= get_edge_psi('qflux_5_300', lonin=lonin, thresh=0., intdown=intdown)
vars_10_300 = get_edge_psi('qflux_10_300', lonin=lonin, thresh=0., intdown=intdown)
vars_15_300 = get_edge_psi('qflux_15_300', lonin=lonin, thresh=0., intdown=intdown)
vars_20_300 = get_edge_psi('qflux_20_300', lonin=lonin, thresh=0., intdown=intdown)
vars_25_300 = get_edge_psi('qflux_25_300', lonin=lonin, thresh=0., intdown=intdown)#

vars_90 = get_edge_psi('ss_90.000', thresh=0., intdown=intdown)
vars_91 = get_edge_psi('ss_91.000', thresh=0., intdown=intdown)
vars_92 = get_edge_psi('ss_92.000', thresh=0., intdown=intdown)
vars_93 = get_edge_psi('ss_93.000', thresh=0., intdown=intdown)
vars_94 = get_edge_psi('ss_94.000', thresh=0., intdown=intdown)
vars_95 = get_edge_psi('ss_95.000', thresh=0., intdown=intdown)
vars_100 = get_edge_psi('ss_100.000', thresh=0., intdown=intdown)
vars_105 = get_edge_psi('ss_105.000', thresh=0., intdown=intdown)


plt.plot(vars_0_100[0], -1.*vars_0_100[1], 'xk', mew=2)
plt.plot(vars_5_100[0], -1.*vars_5_100[1], 'xb', mew=2)
plt.plot(vars_10_100[0], -1.*vars_10_100[1], 'xg', mew=2)
plt.plot(vars_15_100[0], -1.*vars_15_100[1], 'xr', mew=2)
plt.plot(vars_20_100[0], -1.*vars_20_100[1], 'xy', mew=2)
plt.plot(vars_25_100[0], -1.*vars_25_100[1], 'xm', mew=2)

plt.plot(vars_0_200[0], -1.*vars_0_200[1], '+k', mew=2)
plt.plot(vars_5_200[0], -1.*vars_5_200[1], '+b', mew=2)
plt.plot(vars_10_200[0], -1.*vars_10_200[1], '+g', mew=2)
plt.plot(vars_15_200[0], -1.*vars_15_200[1], '+r', mew=2)
plt.plot(vars_20_200[0], -1.*vars_20_200[1], '+y', mew=2)
plt.plot(vars_25_200[0], -1.*vars_25_200[1], '+m', mew=2)

plt.plot(vars_0_300[0], -1.*vars_0_300[1], 'ok', mew=2)
plt.plot(vars_5_300[0], -1.*vars_5_300[1], 'ob', mew=2)
plt.plot(vars_10_300[0], -1.*vars_10_300[1], 'og', mew=2)
plt.plot(vars_15_300[0], -1.*vars_15_300[1], 'or', mew=2)
plt.plot(vars_20_300[0], -1.*vars_20_300[1], 'oy', mew=2)
plt.plot(vars_25_300[0], -1.*vars_25_300[1], 'om', mew=2)

plt.plot(vars_90[0], -1.*vars_90[1], 'sk', mew=2)
plt.plot(vars_91[0], -1.*vars_91[1], 'sk', mew=2)
plt.plot(vars_92[0], -1.*vars_92[1], 'sk', mew=2)
plt.plot(vars_93[0], -1.*vars_93[1], 'sk', mew=2)
plt.plot(vars_94[0], -1.*vars_94[1], 'sk', mew=2)
plt.plot(vars_95[0], -1.*vars_95[1], 'sk', mew=2)
plt.plot(vars_100[0], -1.*vars_100[1], 'sk', mew=2)
plt.plot(vars_105[0], -1.*vars_105[1], 'sk', mew=2)



plt.xlabel('Cell edge')
plt.ylabel('Max 500 hPa Mass Streamfunction')
plt.grid(True,linestyle=':')

plt.xscale('log')
plt.yscale('log')
plt.minorticks_off()
#plt.xlim(0.625,40)
#plt.ylim(100,600)
#plt.xticks([1.25,2.5,5,10,20,40])
#plt.yticks([75,150,300,600])
plt.xlim(1,30)
plt.ylim(100,400)
plt.xticks([1,5,10,20])
plt.yticks([100,200,300,400])
ax=plt.gca()
ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
ax.get_yaxis().set_major_formatter(tk.ScalarFormatter())
plt.tight_layout()

plt.savefig('/scratch/rg419/plots/steady_state_runs/regimefig_island.pdf', format='pdf')
plt.close()
            
            