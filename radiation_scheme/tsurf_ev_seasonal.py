# Plot the evolution of surface temperature as a function of time

import numpy as np
import os
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import cell_area

def tsurf_ev(run, months):
    
    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/atmos_monthly.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    #read data into xarray 
    data = xr.open_mfdataset( names, decode_times=False)
    data.coords['year'] = data.time//360 + 1        
    data_yr = data.groupby('year').mean(('time'))
    
    area = cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')
    area_xr = xr.DataArray(area, [('lat', data.lat ), ('lon', data.lon)])
        
    #take area mean of t_surf
    t_surf_av = (data_yr.t_surf*area_xr).sum(('lat','lon'))/area_xr.sum(('lat','lon'))

    t_surf_av.plot()
    #plt.xlabel('Pentad')
    plt.ylabel('Global mean temperature, K')
    plotname = '/scratch/rg419/plots/radiation_scheme/t_surf_ev_'+ run +'.png'
    plt.savefig(plotname)
    plt.close()
    
    return t_surf_av



if __name__ == "__main__":


	#t_surf_av = tsurf_ev('ap_2', [13,481])
	#t_surf_av = tsurf_ev('gvr_rrtm', [1,481])
	t_surf_av = tsurf_ev('gvr_geen', [1,481])