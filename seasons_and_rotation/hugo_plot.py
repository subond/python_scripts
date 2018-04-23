# Evaluate the mass streamfunction as a function of time, look at progression of peak throughout year

from physics import mass_streamfunction
from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

L = 2.500e6
cp =  287.04/(2./7.)
g = 9.8
stefan = 5.6734e-8

def hugo_plot(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1.):
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    omega_max_loc = data.lat[np.argmin(data.omega[:,27,:,:].mean('lon').values, axis=1)]
    tsurf_max_loc = data.lat[np.argmax(data.t_surf.mean('lon').values, axis=1)]
    
    plt.plot(tsurf_max_loc, omega_max_loc, 'x')
    plt.plot([-50,50],[-50,50],'k--')
    plt.xlim(-50,50)
    plt.ylim(-50,50)
    plt.xlabel('Latitude of max SST')
    plt.ylabel('Latitude of max w')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/hugo_'+run+'.png')
    plt.show()

hugo_plot('ap_2_8core_720', [121,409])

#hugo_plot('ap_1', [121,157])
#hugo_plot('ap_30', [121,409])
#hugo_plot('aquaplanet_10m', [121,481])
#hugo_plot('sn_3.000', [121,481], period_fac=3.)

