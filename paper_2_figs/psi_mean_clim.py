"""
Plot the climatology of the global mean psi to compare with equilibration experiments
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
from hadley_cell import mass_streamfunction
from data_handling_updates import cell_area


def psi_mean_clim(run):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    
    area = cell_area(42, '/scratch/rg419/Isca/')
    area_xr = xr.DataArray(area, [('lat', data.lat ), ('lon', data.lon)])
    area_xr = area_xr.mean('lon')
    psi_mean = ((psi*area_xr).sum(('lat'))/area_xr.sum(('lat'))).mean('pfull')
    psi_mean = psi_mean*-1.

    return psi_mean

# Set plotting directory
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


psi_p25 = psi_mean_clim('sn_0.250')
psi_p5 = psi_mean_clim('sn_0.500')
psi_1 = psi_mean_clim('sn_1.000')
psi_2 = psi_mean_clim('sn_2.000')
psi_4 = psi_mean_clim('sn_4.000')

plt.plot(psi_p25.xofyear*4./5., psi_p25)
plt.plot(psi_p5.xofyear*2., psi_p5)
plt.plot(psi_1.xofyear, psi_1)
plt.plot(psi_2.xofyear/2., psi_2)
plt.plot(psi_4.xofyear/4., psi_4)

plt.grid(True,linestyle=':')
plt.savefig(plot_dir + 'psi_mean_sn.pdf', format='pdf')
plt.close()


psi_2p5 = psi_mean_clim('mld_2.5')
psi_5 = psi_mean_clim('mld_5')
psi_10 = psi_mean_clim('sn_1.000')
psi_15 = psi_mean_clim('mld_15')
psi_20 = psi_mean_clim('mld_20')

plt.plot(psi_2p5.xofyear, psi_2p5)
plt.plot(psi_5.xofyear, psi_5)
plt.plot(psi_10.xofyear, psi_10)
plt.plot(psi_15.xofyear, psi_15)
plt.plot(psi_20.xofyear, psi_20)

plt.grid(True,linestyle=':')
plt.savefig(plot_dir + 'psi_mean_mld.pdf', format='pdf')
plt.close()