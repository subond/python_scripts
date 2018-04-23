# Load in data for steady state simulations and save average of years 11-15 in climatology folder. 

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import cell_area
from physics import mass_streamfunction
import scipy.interpolate as spint

def get_edge_psi_ss(run, lonin=[-1.,361.], sanity_check=False, lev=500., thresh=120., intdown=False):
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    # Calculate mass streamfunction
    psi_coarse = mass_streamfunction(data, a=6376.0e3, dp_in=50., lons=lons, intdown=intdown)
    psi_coarse /= 1.e9
    
    # Interpolate between 35S and 35N (Should help to avoid mistakes with Ferrel cell etc)
    f = spint.interp1d(psi_coarse.lat, psi_coarse.sel(pfull=lev), axis=0, fill_value='extrapolate')
    lats = np.arange(-35, 35.1, 0.5)
    psi = f(lats)
    psi = xr.DataArray(psi, coords=[lats], dims=['lat'])
    
    psi_mask = np.ones(psi.values.shape)
    
    if thresh == 0.:
        psi_mask[psi.values <= thresh] = np.nan
    else:
        psi_mask[np.abs(psi).values <= thresh] = np.nan
    
    psi_mask = xr.DataArray( psi_mask, psi.coords)
    psi_masked = psi * psi_mask
    
    if thresh ==0.:
        psi_red = psi
    else:
        psi_red = psi_masked
    
    psi_mask_red = np.nan_to_num(psi_mask)
    
    psi_max = np.abs(psi_red).max('lat')
    psi_max_loc = psi_red.lat.values[np.abs(psi_red).argmax('lat').values]
    
    # Create list of which times the max of abs(psi_rad) is positive vs negative
    ispos = psi_red.values[np.abs(psi_red).argmax('lat').values] >= 0.
    
    # Use this first to get the sign of psi_max right
    if not ispos: psi_max = -1.*psi_max
    
    # Locate places where the psi mask changes sign    
    edges_are_at = psi_red.lat[np.where(psi_mask_red[:-1] != psi_mask_red[1:])[0]]
    try:
        # Use a try statement to avoid errors where edge is poorly defined
        if ispos:
            # If the overturning is positive, look for the first edge to the south of the max
            edge_loc = np.max(edges_are_at[edges_are_at < psi_max_loc])
        else:
            # If the overturning is negative, look for the first edge to the north of the max
            edge_loc = np.min(edges_are_at[edges_are_at > psi_max_loc])
    except:
        # If we can't find an edge, set the value to nan 
        edge_loc = np.nan
    
    if sanity_check:
        psi_coarse.plot.contourf(x='lat', y='pfull', yincrease=False, levels=np.arange(-500.,501.,50.))
        plt.plot(edge_loc, 500., 'kx', mew=2)
        plt.plot(psi_max_loc, 500., 'kx', mew=2)
        plt.show()

    return edge_loc, psi_max

if __name__ == "__main__":
    get_edge_psi_ss('ss_90.000', sanity_check=True, thresh=0.)
    
    
    