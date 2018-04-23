# Find the peak in moist static energy for a given data file

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import scipy.interpolate as spint


def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons
    
def peak_mse(data, lonin=[-1.,361.], lat_bound=45.):        
    
    lons = pick_lons(data, lonin)
    
    g = 9.8
    cp = 287.04/2*7
    L = 2.500e6
    
    mse = (data.temp*cp + data.sphum*L + data.height*g).sel(lon=lons).mean('lon').sel(pfull=900.)
    
    f = spint.interp1d(mse.lat, mse, axis=-1, fill_value='extrapolate', kind='quadratic')
    lats_new = np.arange(-lat_bound, lat_bound+0.1, 0.1)
    mse = f(lats_new)
    
    # Determine whether p_new is 2d or 1d and create DataArray
    if 'xofyear' in data.coords:
        mse = xr.DataArray(mse, coords=[data.xofyear.values, lats_new], dims=['xofyear', 'lat'])
    else:
        mse = xr.DataArray(mse, coords=[lats_new], dims=['lat'])
    
    mse = mse/1000.
    mse_max = mse.max('lat')
    mse_max_loc = mse.lat.values[mse.argmax('lat').values]
    mse_max_loc = xr.DataArray(mse_max_loc, coords=mse_max.coords, dims=mse_max.dims)
    
    if 'xofyear' in data.coords:
        output_dict = {'mse': (['xofyear', 'lat'], mse),
                   'mse_max': (['xofyear'], mse_max),
                   'mse_max_loc': (['xofyear'], mse_max_loc)}
    
        coord_dict = {'xofyear': ('xofyear', mse.xofyear),
                      'lat': ('lat', mse.lat)}
    else:
        output_dict = {'mse': (['lat'], mse),
                   'mse_max': (mse_max),
                   'mse_max_loc': (mse_max_loc)}
    
        coord_dict = {'lat': ('lat', mse.lat)}
                      
    mse_data = xr.Dataset(output_dict, coords=coord_dict)

    return mse_data
    

if __name__ == "__main__":
    # Sanity check
    lonin = [-1.,361.]
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')
    data = peak_mse(data, lonin=lonin)
        
    if 'xofyear' in data.coords:
        data.mse.plot.contourf(x='xofyear', y='lat')
        data.mse_max_loc.plot.line('k')
    else:
        data.mse.plot.line('k')
        plt.plot(data.mse_max_loc, data.mse_max, 'xr')
    plt.show()
    
    plt.plot(data.mse_max_loc, data.mse_max, 'xk')
    plt.show()
