# Plot timeseries of surface fluxes and SST at 15N

from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

L = 2.500e6
cp =  287.04/(2./7.)
g = 9.8
stefan = 5.6734e-8

def surf_ts(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1.):
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    plt.plot(data.xofyear,data.flux_t[:,37,:].mean('lon'))
    plt.xlabel('Pentad')
    plt.ylabel('Sensible heat flux, W/m2 (15N)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_flux_t_ts.png')
    plt.close()

    plt.plot(data.xofyear,data.flux_lhe[:,37,:].mean('lon'))
    plt.xlabel('Pentad')
    plt.ylabel('Latent heat flux, W/m2 (15N)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_flux_lhe_ts.png')
    plt.close()
    
    plt.plot(data.xofyear,data.flux_sw[:,37,:].mean('lon'))
    plt.xlabel('Pentad')
    plt.ylabel('SW heat flux, W/m2 (15N)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_flux_sw_ts.png')
    plt.close()
    
    bsurf = stefan*np.power(data.t_surf,4)
    
    plt.plot(data.xofyear,(bsurf-data.flux_lw)[:,37,:].mean('lon'))
    plt.xlabel('Pentad')
    plt.ylabel('LW heat flux, W/m2 (15N)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_flux_lw_ts.png')
    plt.close()
    
    plt.plot(data.xofyear,(bsurf-data.flux_lw + data.flux_lhe)[:,37,:].mean('lon'))
    plt.xlabel('Pentad')
    plt.ylabel('LW heat flux, W/m2 (15N)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_flux_lwlhe_ts.png')
    plt.close()
    
    plt.plot(data.xofyear,data.t_surf[:,37,:].mean('lon'))
    plt.xlabel('Pentad')
    plt.ylabel('SST, K (15N)')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_SST15N_ts.png')
    plt.close()
    
    print 'done'
    

surf_ts('aquaplanet_10m', [121,481])

#hm_vars('aquaplanet_2m', [121,481], filename='atmos_daily')

#hm_vars('sn_3.000', [289,469], period_fac=3.)