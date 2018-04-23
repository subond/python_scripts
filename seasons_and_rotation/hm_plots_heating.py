# Evaluate vertically integrated heating terms and plot as fn of time and lat

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
    dp=xr.DataArray(data.phalf.diff('phalf').values*100, coords=[('pfull', data.pfull)])
    
    cndht = (data.dt_tg_condensation.mean('lon')*dp).sum('pfull')
    cnvht = (data.dt_tg_convection.mean('lon')*dp).sum('pfull')
    difht = (data.dt_tg_diffusion.mean('lon')*dp).sum('pfull')
    radht = (data.tdt_rad.mean('lon')*dp).sum('pfull')
    totht = cndht + cnvht + radht + difht
    
    cndht.plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both',levels=np.arange(-0.75,0.8,0.05))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('Condensation Heating K/s*kg/m2')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_cndht_hm.png')
    plt.close()
    
    cnvht.plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both',levels=np.arange(0,6.5,0.5))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('Convection Heating K/s*kg/m2')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_cnvht_hm.png')
    plt.close()
    
    radht.plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both',levels=np.arange(-2.,0.2,0.1))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('Radiative Heating K/s*kg/m2')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_radht_hm.png')
    plt.close()
    
    difht.plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both',levels=np.arange(-0.24,0.26,0.02))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('Diffusive Heating K/s*kg/m2')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_difht_hm.png')
    plt.close()
    
    totht.plot.contourf(x='xofyear',y='lat', add_labels=False, extend='both',levels=np.arange(-6.,6.5,1.))
    plt.grid()
    plt.ylim(-45,45)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title('Total Heating K/s*kg/m2')
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+run+'_totht_hm.png')
    plt.close()
    
    print 'done'
    

hm_vars('aquaplanet_10m', [121,481])
hm_vars('sn_3.000', [121,481], period_fac=3.)



