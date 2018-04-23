"""
Calculate precipitation centroid. Calculation based on Frierson and Hwang 2012 (J. Clim.) and Donohoe et al. 2013 (J. Clim.)

"""

import xarray as xr
from data_handling import cell_area
import scipy.interpolate as spint
import numpy as np

def precip_centroid(run, lat_bound=45., lonin=[-1.,361.]):
    
    area = cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')   # Area of grid cells
    
    #Load in data, add area to dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['area'] = (('lat','lon'), area)
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    # Get total precip
    try:
        data['precipitation'] = data.condensation_rain + data.convection_rain
    except:
        data['precipitation'] = data.precipitation
    
    # Select latitudes over which to evaluate precip centroid
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -lat_bound and data.lat[i] <= lat_bound]
    
    # Integrate precip wrt longitude
    precip_area_lats = (data.precipitation.sel(lat=lats) * data.area.sel(lat=lats)).sel(lon=lons).sum('lon').values
    
    # Interpolate precip in latitude    
    f = spint.interp1d(lats, precip_area_lats, axis=1, fill_value='extrapolate')
    lats_new = np.arange(-lat_bound, lat_bound+0.1, 0.1)
    p_new = f(lats_new)
    p_new = xr.DataArray(p_new, coords=[data.xofyear.values, lats_new], dims=['xofyear', 'lat'])
    
    # Calculate cumulative sum of precip with latitude
    p_area_int = p_new.cumsum('lat')
    
    # At each time find the precipitation centroid: the latitude at which half of the area integrated precip lies North/South
    p_cent = np.zeros((len(p_new.xofyear.values),))
    for i in range(1,len(p_new.xofyear.values)+1):
        p_cent[i-1] = p_new.lat[p_area_int.sel(xofyear=i) <= 0.5 * p_area_int.sel(xofyear=i).max('lat')].max('lat').values
        
    p_cent= xr.DataArray(p_cent, coords=[p_new.xofyear.values], dims=['xofyear'])
    
    return p_cent
    

