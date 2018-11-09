"""
Plot precip and mse hm and overturning/ascent for 10m mld run, plus plot spin up progression of these

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_mse_plot
from pylab import rcParams
from hadley_cell import mass_streamfunction
from data_handling_updates import make_sym

plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
rcParams['figure.figsize'] = 5.5, 3.5
rcParams['font.size'] = 14


fig, (ax1) = plt.subplots(1, sharex=True)

def precip_psi_plot(run, ax, label=''):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
        
    f1 = precip_mse_plot(data, ax, plot_type='precip', precip_contour=None)
    
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    
    ax.text(-10, 60, label)
    
    figname='precip_mse_' + run + '.pdf'
    
    return f1, figname

#subsolar_point = -23.439*np.cos(np.arange(0.,361.,1.)*np.pi/180.)

#f1 = precip_psi_plot('zs_10m', ax=ax1)
#f1, figname = precip_psi_plot('zs_10m_qflux', ax=ax1)
f1, figname = precip_psi_plot('zs_10m_qflux_bs', ax=ax1)

ax1.set_xticks([12,24,36,48,60,72])
ax1.set_xlabel('Pentad')

ax1.set_xlim([1,72])

plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.05)

cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.2, pad=0.2, aspect=40)
cb1.set_label('Precipitation, mm/day')

plt.savefig(plot_dir+figname, format='pdf')
plt.close()        
