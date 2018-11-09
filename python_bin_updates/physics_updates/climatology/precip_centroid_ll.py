"""
Calculate precipitation centroid. Calculation based on Frierson and Hwang 2012 (J. Clim.) and Donohoe et al. 2013 (J. Clim.)
Modified 1/05/2018 to take multiple years of input data to allow error calculations

"""

import xarray as xr
from data_handling_updates import cell_area, gradients as gr
import scipy.interpolate as spint
import numpy as np
import matplotlib.pyplot as plt


    
def precip_centroid_ll(data, lat_bound=45., res=0.01):
    '''Inputs: data = input DataSet
               lat_bound = lat range to integrate over
               lonin = longitude range to use'''
    
    if not 'area' in data.data_vars:
        area = cell_area(42, '/scratch/rg419/Isca/')   # Area of grid cells
        # Add area to dataset
        data['area'] = (('lat','lon'), area)
        
    # Get total precip
    try:
        data['precipitation'] = data.condensation_rain + data.convection_rain
    except:
        data['precipitation'] = data.precipitation
    
    # Select latitudes over which to evaluate precip centroid
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -lat_bound and data.lat[i] <= lat_bound]

    # area weight precip
    precip_area_lats = (data.precipitation.sel(lat=lats) * data.area.sel(lat=lats)).values
    
    # Interpolate precip in latitude    
    lat_ax = data.precipitation.get_axis_num('lat')
    f = spint.interp1d(lats, precip_area_lats, axis=lat_ax, fill_value='extrapolate', kind='quadratic')
    lats_new = np.arange(-lat_bound, lat_bound+res, res)
    p_new = f(lats_new)
    
    # Determine the time dimensions of p_new and create DataArray
    if 'year_no' in data.coords:
        p_new = xr.DataArray(p_new, coords=[data.year_no.values, data.xofyear.values, lats_new, data.lon.values], dims=['year_no', 'xofyear', 'lat', 'lon'])
    elif 'xofyear' in data.coords:
        p_new = xr.DataArray(p_new, coords=[data.xofyear.values, lats_new, data.lon.values], dims=['xofyear', 'lat', 'lon'])
    elif 'time' in data.coords:
        p_new = xr.DataArray(p_new, coords=[data.time.values, lats_new, data.lon.values], dims=['time', 'lat', 'lon'])
    else:
        p_new = xr.DataArray(p_new, coords=[lats_new, data.lon.values], dims=['lat','lon'])
        
    #try:
    #     p_new = xr.DataArray(p_new, coords=[data.xofyear.values, lats_new], dims=['xofyear', 'lat'])
    #except:
    #    p_new = xr.DataArray(p_new, coords=[lats_new], dims=['lat'])
            
    # Calculate cumulative sum of precip with latitude
    p_area_int = p_new.cumsum('lat')
    
    # At each time find the precipitation centroid: the latitude at which half of the area integrated precip lies North/South
    # Loop over time dimensions, if these exist
            
    if 'xofyear' in data.coords:
        p_cent = np.zeros((len(p_new.xofyear.values),len(p_new.lon.values),))
        j=0
        for i in p_new.xofyear.values:
            l=0
            for k in p_new.lon.values:
                #p_area_int.sel(xofyear=i).sel(lon=k).plot.line()
                p_cent[j,l] = p_new.lat[p_area_int.sel(xofyear=i).sel(lon=k) <= 0.5 * p_area_int.sel(xofyear=i).sel(lon=k).max('lat')].max('lat').values
                l=l+1
            j=j+1
        
        p_cent= xr.DataArray(p_cent, coords=[p_new.xofyear.values, p_new.lon.values], dims=['xofyear','lon'])
    
    elif 'time' in data.coords:
        p_cent = np.zeros((len(p_new.time.values),len(p_new.lon.values),))
        j=0
        for i in p_new.time.values:
            l=0
            for k in p_new.lon.values:
                p_cent[j,l] = p_new.lat[p_area_int.sel(time=i).sel(lon=k) <= 0.5 * p_area_int.sel(time=i).sel(lon=k).max('lat')].max('lat').values
                l=l+1
            j=j+1
        
        p_cent= xr.DataArray(p_cent, coords=[p_new.time.values, p_new.lon.values], dims=['time','lon'])
    
    else:
        #print(p_area_int.max('lat'))
        p_cent = np.zeros((len(p_new.lon.values),))
        l=0
        for k in p_new.lon.values:
            p_cent[l] = p_new.lat[p_area_int.sel(lon=k) <= 0.5 * p_area_int.sel(lon=k).max('lat')].max('lat').values
            l=l+1
        
        p_cent= xr.DataArray(p_cent, coords=[p_new.lon.values], dims=['lon'])
        
        
    #try:
    #    p_cent = np.zeros((len(p_new.xofyear.values),))
    #    for i in range(1,len(p_new.xofyear.values)+1):
    #        p_cent[i-1] = p_new.lat[p_area_int.sel(xofyear=i) <= 0.5 * p_area_int.sel(xofyear=i).max('lat')].max('lat').values
        
    #    p_cent= xr.DataArray(p_cent, coords=[p_new.xofyear.values], dims=['xofyear'])
    #except:
    #    p_cent = p_new.lat[p_area_int <= 0.5 * p_area_int.max('lat')].max('lat').values
    
    data['p_cent'] = p_cent
    
    return data
    

if __name__ == "__main__":
    # Sanity check
    import sh
    #data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/qflux_ss/qflux_10_200.nc')
    data_dir = '/disca/restore/gv2scratch/rg419/Data_moist/climatologies/qflux_ss/'
    runs = ['qflux_0_100.nc', 'qflux_0_200.nc', 'qflux_0_300.nc',
            'qflux_5_100.nc', 'qflux_5_200.nc', 'qflux_5_300.nc',
            'qflux_10_100.nc', 'qflux_10_200.nc', 'qflux_10_300.nc',
            'qflux_15_100.nc', 'qflux_15_200.nc', 'qflux_15_300.nc',
            'qflux_20_100.nc', 'qflux_20_200.nc', 'qflux_20_300.nc',
            'qflux_25_100.nc', 'qflux_25_200.nc', 'qflux_25_300.nc']
    
    plot_dir = '/scratch/rg419/plots/qflux_ss/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    for run in runs[15:18]:
        data = xr.open_dataset(data_dir + run)
        data = precip_centroid_ll(data)
        data.p_cent.plot.line()
        
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    plt.savefig(plot_dir+'precip_centroids_25.pdf', format='pdf')
    plt.close()

