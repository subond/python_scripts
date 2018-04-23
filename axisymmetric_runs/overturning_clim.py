"""
PLot a climatology of the mass streamfunction in the axisymmetric prescribed SST run
"""

import numpy as np
import xarray as xr
from physics import mass_streamfunction
from pylab import rcParams
import sh
import matplotlib.pyplot as plt

    
rcParams['figure.figsize'] = 7, 7
rcParams['font.size'] = 15
#rcParams['text.usetex'] = True

data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/zs_sst.nc')

plot_dir = '/scratch/rg419/plots/axisymmetric_runs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
psi = mass_streamfunction(data, a=6376.0e3, dp_in = 50.)
psi /= 1.e9
    
#levels = np.arange(-2.,2.1,0.25)
    
f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
               
plt.set_cmap('RdBu_r')

ax_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]

levels=np.arange(-300.,301.,60.)

for i in range(12):
    
    ax = ax_list[i]
    
    f1 = psi[:,i,:].plot.contourf(x='lat', y='pfull', yincrease=False, levels=levels, ax=ax, extend = 'both', add_labels=False, add_colorbar=False)
    ax.set_xticks(range(-30,31,60))
    ax.set_xlim(-60.,60.)
    ax.grid(True,linestyle=':')


plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.05, hspace=0.25, wspace=0.15)

cb1=f.colorbar(f1, ax=(ax_list), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)

figname = 'psi_monthly.pdf'

plt.savefig(plot_dir + figname, format='pdf')
plt.close()




# Also produce this plot for sn_1.000 for comparison (eddy permitted equivalent)

data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')
data.coords['month'] =  (data.xofyear -1)//6 + 1
data = data.groupby('month').mean(('xofyear'))
psi_ep = mass_streamfunction(data, a=6376.0e3, dp_in = 50.)
psi_ep /= 1.e9

f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
ax_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]

for i in range(12):
    
    ax = ax_list[i]
    
    f1 = psi_ep[:,i,:].plot.contourf(x='lat', y='pfull', yincrease=False, levels=levels, ax=ax, extend = 'both', add_labels=False, add_colorbar=False)
    ax.set_xticks(range(-30,31,60))
    ax.set_xlim(-60.,60.)
    ax.grid(True,linestyle=':')


plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.05, hspace=0.25, wspace=0.15)

cb1=f.colorbar(f1, ax=(ax_list), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)

figname = 'psi_eddies_monthly.pdf'

plt.savefig(plot_dir + figname, format='pdf')
plt.close()
