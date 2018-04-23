#Functions for loading and time averaging netcdf data ready to plot

from netCDF4 import Dataset
import pandas as pd
import xarray as xr


# Load netcdf data for a given run for desired year
def load_year_xr(run_fol, year, pinterp=False, era=True):
    if year == 1:
        #spin up directory
        names = ['/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/' + run_fol + '/run%d/atmos_daily.nc' % m for m in range(1, 13)]
    elif pinterp == True and era ==True:
        #Data directory, pressure interpolated files:
        names = ['/scratch/rg419/Data_moist/' + run_fol + '/run%d/plev_daily.nc' % m for m in range( year*12-11, year*12+1)  ]
    elif pinterp == True and era == False:
        names = ['/scratch/rg419/Data_moist/' + run_fol + '/run%d/pmodel_daily.nc' % m for m in range( year*12-11, year*12+1)  ]
    else:
        #Data directory:
        names = ['/scratch/rg419/Data_moist/' + run_fol + '/run%d/atmos_daily.nc' % m for m in range( year*12-11, year*12+1)  ]

    #read data into xarray 
    rundata = xr.open_mfdataset( names,
      decode_times=False,  # no calendar so tell netcdf lib
      # choose how data will be broken down into manageable chunks.
       chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    return rundata


# Load netcdf data for a given run for desired years
def load_syear_xr(run_fol, year, pinterp=False):

    if pinterp == True:
        #Data directory, pressure interpolated files:
        names = ['/scratch/rg419/Data_moist/' + run_fol + '/run%d/plev_daily.nc' % m for m in range(12*year+3,12*year+15)  ]
    else:
        #Data directory:
        names = ['/scratch/rg419/Data_moist/' + run_fol + '/run%d/atmos_daily.nc' % m for m in range(12*year+3,12*year+15)  ]

    #read data into xarray 
    rundata = xr.open_mfdataset( names,
          decode_times=False,  # no calendar so tell netcdf lib
          # choose how data will be broken down into manageable chunks.
	  chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})

    return rundata




# Load netcdf data for a given season for desired year.  NB. years defined to run from March - February + this script won't handle spin-up data.
def load_season_xr(run_fol, season, year, pinterp=False):
    season_dic = {'mam':range(12*year+3,12*year+6),
                  'jja':range(12*year+6,12*year+9),
                  'son':range(12*year+9,12*year+12),
                  'djf':range(12*year+12,12*year+15)}

    if pinterp == True:
        #Data directory, pressure interpolated files:
        names = ['/scratch/rg419/Data_moist/' + run_fol + '/run%d/plev_daily.nc' % m for m in season_dic[season]  ]
    else:
        #Data directory:
        names = ['/scratch/rg419/Data_moist/' + run_fol + '/run%d/atmos_daily.nc' % m for m in season_dic[season]  ]

    #read data into xarray 
    rundata = xr.open_mfdataset( names,
          decode_times=False,  # no calendar so tell netcdf lib
          # choose how data will be broken down into manageable chunks.
	  chunks={'time': 30, 'lon': 128//4, 'lat': 64//2})
    return rundata



