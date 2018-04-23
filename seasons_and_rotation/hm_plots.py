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

def hm_vars(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1.):
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    psi = mass_streamfunction(data, a=6376.0e3)/1e9
    mse = (cp*data.temp + L*data.sphum + g*data.height)/1000.
    
    mse_max_loc = data.lat[np.argmax(mse[:,38,:,:].mean('lon').values, axis=1)]
    omega_max_loc = data.lat[np.argmin(data.omega[:,27,:,:].mean('lon').values, axis=1)]
    tsurf_max_loc = data.lat[np.argmax(data.t_surf.mean('lon').values, axis=1)]
    
    oei = (np.diff([1 if i>0 else 0 for i in omega_max_loc]).tolist()).index(1)+1
    mei = (np.diff([1 if i>0 else 0 for i in mse_max_loc]).tolist()).index(1)+1
    tei = (np.diff([1 if i>0 else 0 for i in tsurf_max_loc]).tolist()).index(1)+1
    
    data.t_surf.mean('lon').plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(250.,311.,5.))
    plt.grid()
    plt.plot(data.xofyear,tsurf_max_loc)
    plt.plot(data.xofyear, np.amax(tsurf_max_loc) * np.sin((data.xofyear-tei)*5.*np.pi/180./period_fac) )
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('SST, K')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_sst_hm.png')
    plt.close()
    
    mse[:,38,:,:].mean('lon').plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(255,361,5))
    plt.grid()
    plt.plot(data.xofyear,mse_max_loc)
    plt.plot(data.xofyear, np.amax(mse_max_loc) * np.sin((data.xofyear-mei)*5.*np.pi/180./period_fac) )
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('MSE, kJ/kg (921 hPa)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_mse_hm.png')
    plt.close()
    
    data.omega[:,27,:,:].mean('lon').plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(-0.1,0.11,0.01))
    plt.grid()
    plt.plot(data.xofyear,omega_max_loc,'r')
    plt.plot(data.xofyear, np.amax(omega_max_loc) * np.sin((data.xofyear-oei)*5.*np.pi/180./period_fac) ,'g')
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('w, Pa/s (501 hPa)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_w_hm.png')
    plt.close()
    
    psi[:,:,27].plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(-450.,451.,50.))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('Mass streamfunction, 10^9 kg/s (501 hPa)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_psi_hm.png')
    plt.close()
    
    print 'done'
    

#hm_vars('aquaplanet_10m', [121,337])
#hm_vars('sn_3.000', [121,481], period_fac=3.)



def hm_adv(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1.):
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    data.sphum_v[:,38,:,:].mean('lon').plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(-0.15,0.15,0.01))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('vq, kg/kg*m/s (921 hPa)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_vq_hm.png')
    plt.close()
    
    data.vcomp_temp[:,38,:,:].mean('lon').plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(-2000,2001,200))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('vT, Km/s (921 hPa)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_vt_hm.png')
    plt.close()
    
    data.vcomp[:,38,:,:].mean('lon').plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(-9,9.1,0.5))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('v, m/s (921 hPa)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_v_hm.png')
    plt.close()
    
    data.flux_lhe.mean('lon').plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(0,241,20))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('LH, W/m^2 (+ve up)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_lhe_hm.png')
    plt.close()
    
    data.flux_t.mean('lon').plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(-24,24.1,2.))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('SensH, W/m^2 (+ve up)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_sens_hm.png')
    plt.close()
    
    bsurf = stefan*np.power(data.t_surf,4)
    
    (bsurf - data.flux_lw).mean('lon').plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(30,111,5))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('LWsurf, W/m^2 (+ve up)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_lwf_hm.png')
    plt.close()
    
    data.flux_sw.mean('lon').plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both', levels=np.arange(0,351,25))
    plt.grid()
    #plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('SWsurf, W/m^2 (+ve down)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_swf_hm.png')
    plt.close()
    
    print 'done'
    
    

hm_adv('ap_30', [121,409])
hm_vars('ap_30', [121,409])

#hm_adv('aquaplanet_10m', [121,337]) 
#hm_adv('sn_3.000', [121,481], period_fac=3.)
#hm_vars('aquaplanet_2m', [121,481], filename='atmos_daily')
