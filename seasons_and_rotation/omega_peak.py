# Look at progression of peak of omega throughout the year

from physics import mass_streamfunction
from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

a=6376.0e3
omega = 7.2921150e-5

def peak_ascent_psi(run, months, period_fac=1.):
    data = time_means(run, months, filename='atmos_daily', timeav='pentad',period_fac=period_fac)
    psi = mass_streamfunction(data, a=6376.0e3)/1e9
    print 'psi evaluated'
    omega_max_loc = data.lat[np.argmin(data.omega[:,27,:,:].mean('lon').values, axis=1)]
    print 'peak ascent located'

    psi[:,:,36].plot.contourf(x='xofyear', y='lat',add_labels=False)
    plt.plot(data.xofyear, omega_max_loc,'k')
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/wpsi_'+ run +'.png')
    plt.close()
    
#peak_ascent_psi('aquaplanet_2m',[121,481])

#peak_ascent_psi('aquaplanet_10m',[121,481])

#peak_ascent_psi('sn_3.000',[37,73],period_fac=3.)

#peak_ascent_psi('sn_2.000',[25,49],period_fac=2.)


def psi_regime_plots(run, months, period_fac=1.):
    data = time_means(run, months, filename='atmos_daily', timeav='pentad',period_fac=period_fac)
    
    #find index where omega crosses the equator    
    omega_max_loc = data.lat.values[ np.argmin(data.omega[:,27,:,:].mean('lon').values, axis=1) ].tolist()
    oei = (np.diff([1 if i>0 else 0 for i in omega_max_loc]).tolist()).index(1)+1
    
    psi = mass_streamfunction(data.isel(xofyear=[oei-4,oei,oei+4] ), a=6376.0e3)/1e9
    print 'psi evaluated'
    
    u = data.ucomp[[oei-5,oei,oei+5] ,:,:,:].mean('lon')
    a_cos = a*np.cos(data.lat * np.pi/180.)
    m = (u + omega*a_cos)*a_cos
    print 'AM evaluated'
        
    for i in [0,1,2]:
        psi[:,i,:].plot.contourf(x='lat', y='pfull', yincrease=False, add_labels=False, levels=range(-350,351,50))
        m[i,:,:].plot.contour(x='lat', y='pfull', yincrease=False, add_labels=False,colors='k', add_colorbar=False, levels=np.arange(0,np.max(m),omega*a*a/10))
        plt.xlim(-45,45)
        plt.xlabel('Latitude')
        plt.ylabel('Pressure, hPa')
        plt.savefig('/scratch/rg419/plots/seasons_and_rotation/psi_m_'+ run +'_' + str(i) + '.png')
        plt.close()
        print str(i) + ' done'


#psi_regime_plots('aquaplanet_2m',[121,481])
#psi_regime_plots('aquaplanet_10m',[121,481])
#psi_regime_plots('sn_2.000',[25,97],period_fac=2.)

def oei_centred_data(run, months, period_fac=1.):
    data = time_means(run, months, filename='atmos_daily', timeav='pentad',period_fac=period_fac)
    #find index where omega crosses the equator    
    omega_max_loc = data.lat.values[ np.argmin(data.omega[:,27,:,:].mean('lon').values, axis=1) ].tolist()
    oei = (np.diff([1 if i>0 else 0 for i in omega_max_loc]).tolist()).index(1)+1
    
    return omega_max_loc[oei-15:oei+16]
    
oeic_ap2 = oei_centred_data('aquaplanet_2m',[121,481])
oeic_ap10 = oei_centred_data('aquaplanet_10m',[121,481])
oeic_ap10_2 = oei_centred_data('sn_2.000_40',[25,97],period_fac=2.)
oeic_ap10_3 = oei_centred_data('sn_3.000',[37,73],period_fac=3.)

plt.plot(range(-15,16),oeic_ap2,'xk')
plt.plot(range(-15,16),oeic_ap10,'+k')
plt.plot(range(-15,16),oeic_ap10_2,'+r')
plt.plot(range(-15,16),oeic_ap10_3,'+b')
plt.savefig('/scratch/rg419/plots/seasons_and_rotation/quick.png')
    
    
    
    
    