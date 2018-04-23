"""
Get the maximum overturning strength at 500 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from physics import mass_streamfunction, model_constants
import scipy.interpolate as spint

# Define latitudes to interpolate psi onto
lats = np.arange(-35, 35.1, 0.1)


def psi_max_500(run, lonin=[-1.,361.], lev=500., intdown=True):
    """Front end so same function can be used for both climatological and steady state data"""
    
    #Load data
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    # Choose lons to average over
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    # Calculate mass streamfunction
    psi_coarse = mass_streamfunction(data, a=6376.0e3, dp_in=50., lons=lons, intdown=intdown)
    psi_coarse /= 1.e9
    
    # Interpolate between 35S and 35N (Should help to avoid mistakes with Ferrel cell etc)
    f = spint.interp1d(psi_coarse.lat, psi_coarse.sel(pfull=lev), axis=0, fill_value='extrapolate')
    psi = f(lats)
    
    
    # Check if data has a time dimension
    if 'xofyear' in data.coords:
        psi = xr.DataArray(psi, coords=[lats, psi_coarse.xofyear], dims=['lat', 'xofyear'])
        psi_max = xr.DataArray(np.zeros(len(psi.xofyear),), coords=[psi.xofyear], dims=['xofyear'])
        for i in range(len(psi.xofyear)):
            psi_max[i] = max(psi.isel(xofyear=i).max('lat'), psi.isel(xofyear=i).min('lat'), key=abs)
        
    else:
        psi = xr.DataArray(psi, coords=[lats], dims=['lat'])
        psi_max = max(psi.max('lat'), psi.min('lat'), key=abs)

    return psi_max
                            

if __name__ == "__main__":

    psi_max_500('full_qflux', lonin=[60.,150.])
    psi_max_500('ss_95.000')



