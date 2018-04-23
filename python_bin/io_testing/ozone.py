# -*- coding: utf-8 -*-s
import numpy as np
from data_handling import write_to_netcdf as wtn
import xarray as xr


ozone_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/ozone_1990_12.nc'
ozone = xr.open_dataset( ozone_file, decode_times=False)

output_dict={'is_thd':True, 'notime_clim_tser':'clim', 'num_per_year':12, 
             'file_name':'ozone_1990_test.nc', 'var_name':'ozone_1990'}

ozone['ozone_notime'] = (('pfull','lat','lon'), ozone.ozone_1990[0,:,:,:])



wtn.xr_to_nc_file(ozone, 'ozone_1990', output_dict)

