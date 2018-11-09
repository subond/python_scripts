''' 
08/10/2018 Find max or min Walker cell overturning at 500 hPa for various runs
Use similar method to that in hadley_cell/psi_edge_loc.py
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from hadley_cell import walker_cell
from windspharm.xarray import VectorWind
import scipy.interpolate as spint


def walker_strength(data, lonin=[90.,180.], psi_min=True, sanity_check=False):
    '''
    data - a climatology from Isca with dimensions including xofyear
    lonin - longitudes over which to look for min/max values
    psi_min - if true look for minimum of walker streamfunction, otherwise maximum
    sanity_check - produce a plot on screen to check the max/min has been identified correctly
    '''
    
    lons = np.arange(lonin[0], lonin[1], 0.1) # create hi res longitude array to interpolate to
        
    # Create a VectorWind instance, calculate uchi, and use walker_cell to get the streamfunction from this. 
    w = VectorWind(data.ucomp.sel(pfull=np.arange(50.,950.,50.)), data.vcomp.sel(pfull=np.arange(50.,950.,50.)))
    uchi, vchi, upsi, vpsi = w.helmholtz()
    psi_w = walker_cell(uchi, latin=[-20,20], dp_in=-50.)
    # Divide by 1e9 and select 500 hPa level
    psi_w /= 1.e9
    psi_w = psi_w.sel(pfull=500.)
    
    # interpolate to hi res lons
    f = spint.interp1d(psi_w.lon, psi_w, axis=1, fill_value='extrapolate')
    psi_w = f(lons)
    psi_w = xr.DataArray(psi_w, coords=[uchi.xofyear, lons], dims=['xofyear', 'lon'])
    
    # mask to include only pos/neg values
    psi_mask = np.ones(psi_w.values.shape)
    if psi_min:
        psi_mask[psi_w > 0.] = np.nan
    else:
        psi_mask[psi_w < 0.] = np.nan  
    psi_masked = psi_w * psi_mask
    
    # At some times there will be no pos/neg value. Find times where there are lons with non-NaN values
    times_defd = []
    for i in range(0, len(psi_masked.xofyear)):
        if np.any(np.isfinite(psi_masked[i,:])):
            times_defd.append(np.float(psi_masked.xofyear[i]))
            
    # reduce masked array times to include only non-Nan values
    psi_red = psi_masked.sel(xofyear=times_defd)
    
    # Find min/max and its location
    if psi_min:
        walker_mag = psi_red.min('lon')
        walker_mag_loc = psi_red.lon.values[psi_red.argmin('lon').values]
    else:
        walker_mag = psi_red.max('lon')
        walker_mag_loc = psi_red.lon.values[psi_red.argmax('lon').values]
        
    walker_mag_loc = xr.DataArray(walker_mag_loc, coords=[psi_red.xofyear], dims=['xofyear'])
    
    # plot the streamfunction and the location of the maximum as a sanity check
    if sanity_check:
        psi_w.plot.contourf(x='xofyear',y='lon')
        walker_mag_loc.plot()
        #plt.figure(2)
        #plt.plot(walker_mag_loc,walker_mag)
        plt.show()
    
    return walker_mag

if __name__ == "__main__":
    # Use half_shallow as a test case
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/half_shallow.nc')
    walker_strength(data, psi_min=True, lonin=[0.,360.], sanity_check=True)



