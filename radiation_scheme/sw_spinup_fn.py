#Load in the radiative fluxes and look at how these varies with time - is there a signal from the strat wv?

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import cell_area, load_year_xr

def flux_spinup_fn(run_fol,years):
    year = years[0]
    rundata = load_year_xr(run_fol, year)

    area = cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')
    area_xr = xr.DataArray(area, [('lat', rundata.lat ), ('lon', rundata.lon)])
    
    #Initialise arrays to load into
    sw_av =     xr.DataArray(np.zeros((len(years))), [ ('year', years )])
    lw_av =     xr.DataArray(np.zeros((len(years))), [ ('year', years )])

    for year in years:
        print year
        rundata = load_year_xr(run_fol, year)
        sw = rundata.flux_sw.mean(('time'))
        lw = rundata.flux_lw.mean(('time'))

        #take area mean
        sw_in = sw*area_xr
        sw_av[year-years[0]] = sw_in.sum(('lat','lon'))/area_xr.sum(('lat','lon'))
        lw_in = lw*area_xr
        lw_av[year-years[0]] = lw_in.sum(('lat','lon'))/area_xr.sum(('lat','lon'))

    #plot timeseries of these
    plt.figure(1)
    plt.plot(sw_av)
    plt.xlabel('Year')
    plt.ylabel('Mean surface SW flux')
    plt.savefig('/scratch/rg419/plots/sw_spinup.png')
    plt.clf()
    
    plt.figure(2)
    plt.plot(lw_av)
    plt.xlabel('Year')
    plt.ylabel('Mean surface LW flux')
    plt.savefig('/scratch/rg419/plots/lw_spinup.png')
    plt.clf()
    
    return 



if __name__ == "__main__":

        #set run name
        run_fol = 'aquaplanet_10m'

	flux_spinup_fn(run_fol, range(1,41))
