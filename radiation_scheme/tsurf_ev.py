# Plot the evolution of surface temperature as a function of time

import numpy as np
import os
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import cell_area

def tsurf_ev(run, months):
    
    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/atmos_daily.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    #read data into xarray 
    data = xr.open_mfdataset( names, decode_times=False, chunks={'time': 30})
    
    area = cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')
    area_xr = xr.DataArray(area, [('lat', data.lat ), ('lon', data.lon)])
        
    #take area mean of t_surf
    t_surf_av = (data.t_surf*area_xr).sum(('lat','lon'))/area_xr.sum(('lat','lon'))

    t_surf_av.plot()
    plt.xlabel('Day')
    plt.ylabel('Global mean temperature, K')
    #plt.ylim([284,290])
    plotname = '/scratch/rg419/plots/radiation_scheme/t_surf_ev_'+ run +'.png'
    plt.savefig(plotname)
    plt.close()
    
    return t_surf_av



if __name__ == "__main__":

	#t_surf_av = tsurf_ev('vic_bog', [1,61])
	#t_surf_av = tsurf_ev('window_bog', [1,61])
	#t_surf_av = tsurf_ev('basic_bog', [1,61])
	#t_surf_av = tsurf_ev('geen_window', [1,61])
	#t_surf_av = tsurf_ev('geen_nonwindow', [1,61])
	#t_surf_av = tsurf_ev('geen_both', [1,2])
	#t_surf_av = tsurf_ev('geen_log_form', [1,61])
	#t_surf_av = tsurf_ev('geen_log_form_2', [1,28])
	#t_surf_av = tsurf_ev('geen_power_law', [1,2])
	#t_surf_av = tsurf_ev('geen_log_040', [1,61])
	#t_surf_av = tsurf_ev('geen_log_075', [1,61])
	#t_surf_av = tsurf_ev('geen_log_150', [1,61])
	#t_surf_av = tsurf_ev('geen_log_314', [1,61])
	#t_surf_av = tsurf_ev('geen_log_010', [1,121])
	#t_surf_av = tsurf_ev('geen_log_010_init_cond', [1,241])
#	t_surf_av = tsurf_ev('geen_log_040_init_cond', [1,481])
	t_surf_av = tsurf_ev('geen_log_040_init_cond_4co2', [1,450])
	t_surf_av = tsurf_ev('geen_log_040_init_cond_seasons', [1,445])
