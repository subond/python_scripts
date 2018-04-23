#Functions for loading and time averaging netcdf data ready to plot

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap



# Load netcdf data for a given variable for desired months
def multi_month_load(var_name, run_fol, mn_range, pinterp=False, twodim=False):
# Allocate an array depending on whether it's a 2D or 3D variable
    if twodim == True:
        var = np.zeros((len(mn_range),30,64,128))
    elif pinterp == True:
        var = np.zeros((len(mn_range),30,17,64,128))
    else:
        var = np.zeros((len(mn_range),30,40,64,128))
# Loop over the requested months
    i=0
    for mn in mn_range:
# Find data file, either in spin up or Data folder
        if mn < 13:
            nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/' + run_fol + '/run' + str(mn) + '/atmos_daily.nc'
        else:
            if pinterp==True:
                nc_file = '/scratch/rg419/Data_moist/' + run_fol + '/run' + str(mn) + '/plev_daily.nc'
            else:
                nc_file = '/scratch/rg419/Data_moist/' + run_fol + '/run' + str(mn) + '/atmos_daily.nc'
        fh = Dataset(nc_file, mode='r')

        if twodim == True:
            var[i,:,:,:] = fh.variables[var_name][:]
        else:
            var[i,:,:,:,:] = fh.variables[var_name][:]
            var[var>1e20]=np.nan
            #print var[i,10,3,45,:]
        i=i+1
    return var



# Take the average over given individual months
def monthly_mean(var_name, run_fol, mn_range, pinterp=False, twodim=False):
# Use multi month load to load up months
    var = multi_month_load(var_name, run_fol, mn_range, pinterp, twodim)
# Average over days in month
    var = np.mean(var,1)
    return var



# Take the average over multiple months
def intermonthly_mean(var_name, run_fol, mn_range, pinterp=False, twodim=False):
# Load in data using multi month load
    var = multi_month_load(var_name, run_fol, mn_range, pinterp, twodim)
# Average over days in month and months
    var = np.mean(np.mean(var,1),0)
    return var
