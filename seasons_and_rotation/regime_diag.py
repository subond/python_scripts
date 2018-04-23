# Evaluate the mass streamfunction as a function of time, look at progression of peak throughout year

from physics import mass_streamfunction
from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

def regime_diag(run,months,period_fac=1.,rmean=6):
    data = time_means(run, months, filename='atmos_daily', timeav='pentad',period_fac=period_fac)
    psi = mass_streamfunction(data, a=6376.0e3)/1e9
    print 'psi evaluated'

    psi_max = np.argmin( psi.values[0:40,:,:], axis=0)
    
    psi_rate = np.gradient(psi_max[:,27])
    
    def running_mean(x, N):
        cumsum = np.cumsum(np.insert(x, 0, 0)) 
        return (cumsum[N:] - cumsum[:-N]) / N
    
    psi_rate_out = running_mean(psi_rate,rmean)  
    
    return psi_rate_out
    

rate_ap2_1 = regime_diag('aquaplanet_2m',[361,481])
rate_ap10_1 = regime_diag('aquaplanet_10m',[361,481])
rate_ap10_2 = regime_diag('sn_2.000',[25,49],period_fac=2.,rmean=12)
rate_ap10_3 = regime_diag('sn_3.000',[37,73],period_fac=3.,rmean=18)

plt.figure(1)
plt.plot(rate_ap2_1,'kx')
plt.ylim(-1.,1.)
plt.xlabel('Pentad')
plt.ylabel('Rate of change of peak overturning lat')
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/psi_rate_ap2_1.png')
plt.close()

plt.figure(2)
plt.plot(rate_ap10_1,'kx')
plt.ylim(-1.,1.)
plt.xlabel('Pentad')
plt.ylabel('Rate of change of peak overturning lat')
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/psi_rate_ap10_1.png')
plt.close()

plt.figure(3)
plt.plot(rate_ap10_3,'kx')
plt.ylim(-1.,1.)
plt.xlabel('Pentad')
plt.ylabel('Rate of change of peak overturning lat')
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/psi_rate_ap10_3.png')
plt.close()

plt.figure(3)
plt.plot(rate_ap10_2,'kx')
plt.ylim(-1.,1.)
plt.xlabel('Pentad')
plt.ylabel('Rate of change of peak overturning lat')
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/psi_rate_ap10_2.png')
plt.close()

