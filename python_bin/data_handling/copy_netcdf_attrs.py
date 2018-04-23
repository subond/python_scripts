'''Function to copy netcdf attributes to a new files'''

import numpy as np
import xarray as xr
import os
from netcdftime import utime
from netCDF4 import Dataset, date2num


def copy_netcdf_attrs(dsin, dsout, copy_vars = True):
    
    #Copy dimensions
    for dname, the_dim in dsin.dimensions.iteritems():
        print dname, len(the_dim)
        dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

    if copy_vars:
        # Copy variables
        for v_name, varin in dsin.variables.iteritems():
            outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
            print varin.datatype
    
            # Copy variable attributes
            outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    
            outVar[:] = varin[:]
    
    return dsout #returns without closing

