import numpy as np
from data_handling import time_means
import xarray as xr
from stephen import ana, aav, io, sagp
import matplotlib.pyplot as plt


def ape_calc(grid, model_params):
    
    xv, yv = np.meshgrid(grid.lon, grid.lat, indexing='xy')
    grid['lats'] = (('lat','lon'), yv)
        
    ape_sst = 27 * (1 - np.sin(90./60. * grid.lats * np.pi /180.)**2.)
    ape_sst = ape_sst * (np.abs(grid.lats) < 60.)
    ape_sst += 273.15
    ape_sst = np.tile(ape_sst,(12,1,1))
    
    grid['ape_sst'] = (('time','lat','lon'), ape_sst)    
    
    output_dict={'manual_grid_option':False, 'is_thd':False, 'num_years':1., 'time_spacing_days':12, 'file_name':'ape_control_sst.nc', 'var_name':'ape_control_sst'}
    
    io.output_nc_file(grid,'ape_sst', model_params, output_dict)


if __name__ == "__main__":
    
    base_dir= '/scratch/rg419/GFDL_model/GFDLmoistModel/'
    
    resolution_file = base_dir + 'src/extra/python/scripts/gfdl_grid_files/t42.nc'
    
    grid = xr.open_dataset(resolution_file,decode_times=False)
    
    model_params = sagp.model_params_set('/scratch/rg419/GFDL_model/GFDLmoistModel/', delta_t=720., ml_depth=2.)
    
    ape_calc(grid, model_params)


