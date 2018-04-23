# Evaluate the mass streamfunction as a function of time, look at progression of peak throughout year

from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

def omega_max(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1.):
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    omega_max_loc = data.lat[np.argmin(data.omega[:,27,:,:].mean('lon').values, axis=1)]    
    oei = (np.diff([1 if i>0 else 0 for i in omega_max_loc]).tolist()).index(1)+1
    return omega_max_loc, oei

    
oml_ap10, oei_ap10 = omega_max('aquaplanet_10m', [121,481])
oml_ap30, oei_ap30 = omega_max('sn_3.000', [121,481], period_fac=3.)

plt.plot(oml_ap10[oei_ap10-35:oei_ap10+36])
plt.plot(oml_ap30[oei_ap30-35:oei_ap30+36])
plt.ylabel('Latitude')
plt.title('Latitude of peak ascent')
plt.legend(['10m aquaplanet', '30m aquaplanet (3xyear)'],loc='upper left')
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/omega_lat_unscaled.png')
plt.close()


plt.plot(np.arange(1,73),oml_ap10)
plt.plot(np.arange(1.,73,1./3.),oml_ap30)
plt.xlabel('Pentad')
plt.ylabel('Latitude')
plt.title('Latitude of peak ascent')
plt.legend(['10m aquaplanet', '30m aquaplanet (3xyear)'],loc='upper left')
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/omega_lat.png')
plt.close()