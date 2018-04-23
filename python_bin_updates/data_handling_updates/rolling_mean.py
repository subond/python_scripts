import numpy as np
import xarray as xr

def rolling_mean(data, window, tcoord='xofyear'):
    """Perform rolling mean using xarray. NB assumes normal GFDL-MiMA climatology data shape"""
    
    i = window - 1
    
    # Put the last window period in the climatology at the front of the data array
    xofyear_exp = np.zeros([data.shape[0] + i,])
    xofyear_exp[i:] = data.coords[tcoord]
    xofyear_exp[:i] = data.coords[tcoord][-i:]
    
    coords = [('xofyear', xofyear_exp), ('lat', data.lat.values)]
    
    data_exp = xr.DataArray(np.zeros([data.shape[0] + i, data.shape[1]]) , coords=coords)
    data_exp[i:,:] = data
    data_exp[:i,:] = data[-i:,:]
    
    #Take the rolling mean and cut off the NaNs at the front
    data_rm = data_exp.rolling(xofyear=window).mean()
    data_rm = data_rm[i:,:]
    
    return data_rm
    

