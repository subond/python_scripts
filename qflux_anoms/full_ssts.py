# Get SST anoms from full experiment so can base qfluxes on structure

import matplotlib.pyplot as plt
import os
import xarray as xr
import numpy as np
from pylab import rcParams
import sh

rcParams['figure.figsize'] = 10, 6
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

plot_dir = '/scratch/rg419/plots/qflux_anoms/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/full_qflux.nc')
t_anom = data.t_surf - data.t_surf.mean('lon')

t_anom.coords['month'] = (data.xofyear - 1.) //6 +1.
t_anom = t_anom.groupby('month').mean(('xofyear'))
        
for i in range(len(t_anom.month)):
    t_anom.isel(month=i).plot.contourf(x='lon', y='lat', levels=np.arange(-10.,10.,1.), add_labels=False)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.yticks([-90.,-75.,-60.,-45.,-30.,-15.,0.,15.,30.,45.,60.,75.,90.])
    plt.xticks([0.,60.,120.,180.,240.,300.,360.])
    plt.grid(True, linestyle=':')
    figname = 'full_qflux_sst_month' + str(i) + '.pdf'
    plt.savefig(plot_dir + figname)
    plt.close()


p_anom = data.ps - data.ps.mean('lon')

p_anom.coords['month'] = (data.xofyear - 1.) //6 +1.
p_anom = p_anom.groupby('month').mean(('xofyear'))/100.
        
#for i in range(len(t_anom.month)):
#    p_anom.isel(month=i).plot.contourf(x='lon', y='lat', levels=np.arange(-200.,200.,10.), add_labels=False)
#    plt.xlabel('Longitude')
#    plt.ylabel('Latitude')
#    figname = 'full_qflux_ps_month' + str(i) + '.pdf'
#    plt.savefig(plot_dir + figname)
#    plt.close()


(t_anom.isel(month=6) - t_anom.isel(month=0)).plot.contourf(x='lon', y='lat', levels=np.arange(-10.,10.,1.), add_labels=False)
plt.show()