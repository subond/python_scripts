#Functions for loading and time averaging netcdf data ready to plot

from netCDF4 import Dataset
import pandas as pd
import xarray as xr


# Load netcdf data for a given run for desired year
def load_year_xr(run_fol, year, pinterp=False, era=True):
    if pinterp == True and era == True:
        names = ['/scratch/rg419/Data_moist/' + run_fol + '/run%03d/plev_daily.nc' % m for m in range( year*12-11, year*12+1)  ]
    elif pinterp == True and era == False:
        names = ['/scratch/rg419/Data_moist/' + run_fol + '/run%03d/pmodel_daily.nc' % m for m in range( year*12-11, year*12+1)  ]
    else:
        names = ['/scratch/rg419/Data_moist/' + run_fol + '/run%03d/atmos_daily.nc' % m for m in range( year*12-11, year*12+1)  ]
        
    #read data into xarray 
    rundata = xr.open_mfdataset( names,
      decode_times=False,  # no calendar so tell netcdf lib
      # choose how data will be broken down into manageable chunks.
       chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    return rundata


   
if __name__ == "__main__":

    #set run name
        run_fol = 'flat_10m'
        rundata = load_year_xr(run_fol,2)
        print rundata