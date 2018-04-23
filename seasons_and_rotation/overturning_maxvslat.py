# Evaluate the mass streamfunction as a function of time, look at progression of peak throughout year

from physics import mass_streamfunction
from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

def psi_max_lat(run,months,filename='atmos_pentad',period_fac=1.):
    data = time_means(run, months, filename=filename, timeav='pentad',period_fac=period_fac)
    #t_surf_max = data.lat[np.argmax(data.t_surf.mean('lon').values, axis=1)]
    #data.t_surf.mean('lon').plot.contourf(x='xofyear',y='lat')
    omega_max_loc = data.lat[np.argmin(data.omega[36:71,27,:,:].mean('lon').values, axis=1)]
    #plt.ylim(-45,45)
    #plt.plot(t_surf_max)
    
    plt.figure(2)    
    psi = mass_streamfunction(data, a=6376.0e3)/1e9
    print 'psi evaluated'
    #psi_max[0,:] = data.lat[np.argmin( psi.values[:,30:71,27], axis=0)]
    psi_max = np.amin( psi[:,36:71,27], axis=0)
    plt.plot(omega_max_loc[0:11],psi_max[0:11],'x')
    plt.plot(omega_max_loc[12:23],psi_max[12:23],'xr')
    plt.plot(omega_max_loc[24:71],psi_max[24:71],'xg')
    plt.show()


def w_max_lat(run, months,filename='atmos_pentad', period_fac=1.):
    data = time_means(run, months, filename=filename, timeav='pentad',period_fac=period_fac)
    omega_max_loc = data.lat[np.argmin(data.omega[36:71,27,:,:].mean('lon').values, axis=1)]
    omega_max = np.amin(data.omega[36:71,27,:,:].mean('lon').values, axis=1)
    print 'peak ascent located'
    plt.plot(omega_max_loc[0:11],omega_max[0:11],'x')
    plt.plot(omega_max_loc[12:23],omega_max[12:23],'xr')
    plt.plot(omega_max_loc[24:35],omega_max[24:35],'xg')
    plt.show()

psi_max_lat('aquaplanet_10m',[25,301])



