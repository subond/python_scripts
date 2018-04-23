# Evaluate the mass streamfunction as a function of time, look at progression of peak throughout year

from physics import mass_streamfunction
from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

L = 2.500e6
cp =  287.04/(2./7.)
g = 9.8

def hm_vars(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1.):
    
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%d/'+filename+'.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    data = xr.open_mfdataset( names,decode_times=False, chunks={'time': 30})
    
    mse = cp*data.temp + L*data.sphum + g*data.height
    
    mse_max_loc = data.lat[np.argmax(mse[:,38,:,:].mean('lon').values, axis=1)]
    omega_max_loc = data.lat[np.argmin(data.omega[:,27,:,:].mean('lon').values, axis=1)]

    oei = (np.diff([1 if i>0 else 0 for i in omega_max_loc]).tolist()).index(1)+1
    mei = (np.diff([1 if i>0 else 0 for i in mse_max_loc]).tolist()).index(1)+1
    
    plt.plot(data.time,mse_max_loc)
    plt.plot(data.time, -1.*np.amax(mse_max_loc) * np.sin((data.time-mei)*np.pi/180./period_fac) )
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_mse_ts.png')
    plt.close()
    
    plt.plot(data.time,omega_max_loc)
    plt.plot(data.time, -1.*np.amax(omega_max_loc) * np.sin((data.time-oei)*np.pi/180./period_fac) )
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_w_ts.png')
    plt.close()
    
    print 'done'
    

hm_vars('aquaplanet_10m', [241,301])

#hm_vars('aquaplanet_2m', [121,481], filename='atmos_daily')

hm_vars('sn_3.000', [289,469], period_fac=3.)