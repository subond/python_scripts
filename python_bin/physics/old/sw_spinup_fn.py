#Load in the radiative fluxes and look at how these varies with time - is there a signal from the strat wv?

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import cell_area, load_year_xr

def flux_spinup_fn(run, months, filename='atmos_pentad'):

    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    #read data into xarray 
    data = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    
    data.coords['year'] = data.time//360 + 1        
    data_yr = data.groupby('year').mean(('time'))

    area = cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')
    area_xr = xr.DataArray(area, [('lat', data.lat ), ('lon', data.lon)])
    
    #take area mean of fluxes
    sw_av = (data_yr.flux_sw*area_xr).sum(('lat','lon'))/area_xr.sum(('lat','lon'))
    lw_av = (data_yr.flux_lw*area_xr).sum(('lat','lon'))/area_xr.sum(('lat','lon'))
        
    #plot timeseries of these
    plt.figure(1)
    sw_av.plot()
    plt.xlabel('Year')
    plt.ylabel('Mean surface SW flux')
    plt.savefig('/scratch/rg419/plots/spinup/sw_spinup_'+ run +'.png')
    plt.clf()
    
    plt.figure(2)
    lw_av.plot()
    plt.xlabel('Year')
    plt.ylabel('Mean surface LW flux')
    plt.savefig('/scratch/rg419/plots/spinup/lw_spinup_'+ run +'.png')
    plt.clf()
    
    return 



if __name__ == "__main__":


	flux_spinup_fn('ss_90.000', [61,181])
	flux_spinup_fn('ss_95.000', [61,181])
	flux_spinup_fn('ss_100.000', [61,181])
