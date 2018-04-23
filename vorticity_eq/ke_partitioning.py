"""
Look at partitioning between rotational and irrotational KE

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, cell_area
from windspharm.xarray import VectorWind

plot_dir='/scratch/rg419/plots/vorticity_eq/'

def ke_partition(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1.):
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    totp = (data.convection_rain + data.condensation_rain)*86400.

    cell_ar = cell_area(42,'/scratch/rg419/GFDL_model/GFDLmoistModel/')
    cell_ar = xr.DataArray( cell_ar, [('lat', data.lat), ('lon', data.lon )])
    uwnd = data.ucomp
    vwnd = data.vcomp

    w = VectorWind(uwnd, vwnd)

    uchi, vchi, upsi, vpsi = w.helmholtz()

    lats = [i for i in range(len(data.lat)) if data.lat[i] >= 5. and data.lat[i] < 30.]
    lons = [i for i in range(len(data.lon)) if data.lon[i] >= 60. and data.lon[i] < 150.]

    ke_chi = 0.5*(uchi*uchi + vchi*vchi)*cell_ar
    ke_chi_av = ke_chi[:,:,lats,lons].sum('lat').sum('lon')/cell_ar[lats,lons].sum('lat').sum('lon')

    ke_psi = 0.5*(upsi*upsi + vpsi*vpsi)*cell_ar
    ke_psi_av = ke_psi[:,:,lats,lons].sum('lat').sum('lon')/cell_ar[lats,lons].sum('lat').sum('lon')

    totp_av = (totp[:,lats,lons]*cell_ar[lats,lons]).sum('lat').sum('lon')/cell_ar[lats,lons].sum('lat').sum('lon')

    ke_chi_av[:,36].plot()
    ke_psi_av[:,36].plot()
    plt.legend(['KE_chi','KE_psi'])
    plt.xlabel('Pentad')
    plt.ylabel('Kinetic Energy, m2/s2')
    plt.savefig(plot_dir+'KE_'+run+'.png')
    plt.close()

    totp_av.plot()
    plt.xlabel('Pentad')
    plt.ylabel('Precipitation, mm/day')
    plt.savefig(plot_dir+'precip_'+run+'.png')
    plt.close()

ke_partition('ap_1_ed', [121,157])
#ke_partition('era_amipsst', [121,337])

#ke_partition('ap_2_8core_720', [121,409])
#ke_partition('ap_30_st', [193,313])
#ke_partition('ap_1_partraw', [193,313])