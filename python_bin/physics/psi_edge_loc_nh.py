"""
Locate the division between the Hadley cells. Only consider max of SH Hadley cell. 
NB The SS part has not yet been updated here
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from physics import mass_streamfunction, model_constants, gradients as gr

import scipy.interpolate as spint

# Define latitudes to interpolate psi onto
lats = np.arange(-35, 35.1, 0.1)


def get_edge_psi_nh(run, lonin=[-1.,361.], sanity_check=False, lev=500., thresh=0., do_month_av=False, period_fac=1., intdown=True, do_ageostrophic=False, do_geostrophic=False):
    """Front end so same function can be used for both climatological and steady state data"""
    
    #Load data
    if run=='era':
        data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_v_alllevs.nc')
        # Take pentad means
        data.coords['xofyear'] = data.day_of_yr //5 + 1.  
        data = data.groupby('xofyear').mean(('day_of_yr'))
    else:
        data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if do_ageostrophic:
        omega = 7.2921150e-5
        f = 2 * omega * np.sin(data.lat *np.pi/180)
        dphidx = gr.ddx(data.height)
        data = data
        data['vcomp'] = data.vcomp - data.vcomp.mean('lon') - dphidx*9.8/f
        

    if do_geostrophic:
        omega = 7.2921150e-5
        f = 2 * omega * np.sin(data.lat *np.pi/180)
        dphidx = gr.ddx(data.height)
        data['vcomp'] = dphidx*9.8/f 
        
            
    # Choose lons to average over
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    # If data has a time dimension, use time series fn. Otherwise use steady state fn
    if 'xofyear' in data.coords:
        edge_loc, psi_max = get_edge_psi_time(data, lons, sanity_check=sanity_check, 
                             lev=lev, thresh=thresh, intdown=intdown)
                             
        # Take monthly average if wanted, return outputs
        if do_month_av:    
            psi_max.coords['month'] = np.mod( psi_max.xofyear -1, 72.*period_fac) //(6*period_fac) +1    
            psi_max_month = psi_max.groupby('month').mean(('xofyear'))
            edge_loc.coords['month'] = np.mod( edge_loc.xofyear -1, 72.*period_fac) //(6*period_fac) +1    
            edge_loc_month = edge_loc.groupby('month').mean(('xofyear'))
        
            return edge_loc, psi_max, edge_loc_month, psi_max_month
        else:
            return edge_loc, psi_max
            
    else:
        edge_loc, psi_max = get_edge_psi_ss(data, lons, sanity_check=sanity_check, 
                             lev=lev, thresh=thresh, intdown=intdown)
        return edge_loc, psi_max
                             



def get_edge_psi_time(data, lons, sanity_check=False, lev=500., thresh=0., intdown=True):
    """Get cell edge for climatological data"""
    
    # Calculate mass streamfunction
    psi_coarse = mass_streamfunction(data, a=6376.0e3, dp_in=50., lons=lons, intdown=intdown)
    psi_coarse /= 1.e9
    
    # Interpolate between 35S and 35N (Should help to avoid mistakes with Ferrel cell etc)
    f = spint.interp1d(psi_coarse.lat, psi_coarse.sel(pfull=lev), axis=0, fill_value='extrapolate')
    psi = f(lats)
    psi = xr.DataArray(psi, coords=[lats, psi_coarse.xofyear], dims=['lat', 'xofyear'])
    
    # Create a mask that puts NaNs where abs(psi) is below a threshold
    psi_mask = np.ones(psi.values.shape)
    
    if thresh == 0.:
        psi_mask[psi.values <= thresh] = np.nan
    else:
        psi_mask[np.abs(psi).values <= thresh] = np.nan

    psi_mask = xr.DataArray( psi_mask, psi.coords)

    # Use mask to mask psi at 500hPa 
    psi_masked = psi * psi_mask

    # Get a list of times where abs(psi) has values above the threshold
    times_defd = []
    for i in range(0, len(psi_masked.xofyear)):
        if np.any(np.isfinite(psi_masked[:,i])):
            times_defd.append(np.float(psi_masked.xofyear[i]))
    #print times_defd
    
    # Reduce psi to only include times where it goes above the threshold
    if thresh ==0.:
        psi_red = psi.sel(xofyear=times_defd)
    else:
        psi_red = psi_masked.sel(xofyear=times_defd)
        
    psi_mask_red = np.nan_to_num(psi_mask.sel(xofyear=times_defd))
    
    # Up to this point code is same as psi_edge_loc.py. Now take psi_red and look for min instead of taking abs and finding max    
    #psi_max = np.abs(psi_red).max('lat')
    #psi_max_loc = psi_red.lat.values[np.abs(psi_red).argmax('lat').values]
    psi_max = psi_red.min('lat')
    psi_max_loc = psi_red.lat.values[psi_red.argmin('lat').values]
    psi_max_loc = xr.DataArray(psi_max_loc, coords=[psi_red.xofyear], dims=['xofyear'])
    
    
    # _nh Create list of which times the min of psi_rad is positive vs negative
    ispos = []
    for i in range(0, len(psi_max_loc)):
      #  ispos.append(psi_red.values[np.abs(psi_red).argmax('lat').values[i],i] >= 0.)
        ispos.append(psi_red.values[psi_red.argmin('lat').values[i],i] >= 0.)
    
    # Multiply by -1 to give a positive sign for the SH cell
    #psi_max[ispos] = -1.*psi_max
    psi_max = -1.*psi_max
    
    # Now locate cell edge finally!
    edge_loc = np.zeros([len(psi_max_loc),])
    
    for i in range(0, len(psi_max_loc)):
        # Locate places where the psi mask changes sign    
        edges_are_at = psi_red.lat[np.where(psi_mask_red[:-1,i] != psi_mask_red[1:,i])[0]]
        try:
        # Use a try statement to avoid errors where edge is poorly defined
            if ispos[i]:
                # If the overturning is positive, ignore point, as this is from NH cell
                edge_loc[i] = np.nan #np.max(edges_are_at[edges_are_at < psi_max_loc[i]])
            else:
                # If the overturning is negative, look for the first edge to the north of the max
                edge_loc[i] = np.min(edges_are_at[edges_are_at > psi_max_loc[i]])
        except:
        # If we can't find an edge, set the value to nan for this time
            edge_loc[i] = np.nan
        
        
    # Convert to an xarray dataarray for ease of plotting etc.
    edge_loc = xr.DataArray(edge_loc, coords=[psi_red.xofyear], dims=['xofyear'])
    
    # Option to plot up psi as a sanity check that all points are ok
    if sanity_check:
        psi.plot.contourf(levels=np.arange(-500.,510.,50.))
        psi_max_loc.plot()
        edge_loc.plot()
        
        plt.figure(2)
        plt.plot(edge_loc, psi_max, 'x')
        
        plt.show()
    
    return edge_loc, psi_max






def get_edge_psi_ss(data, lons, sanity_check=False, lev=500., thresh=0., intdown=True):
    """Get cell edge for steady state data"""
    
    # Calculate mass streamfunction
    psi_coarse = mass_streamfunction(data, a=6376.0e3, dp_in=50., lons=lons, intdown=intdown)
    psi_coarse /= 1.e9
    
    # Interpolate between 35S and 35N (Should help to avoid mistakes with Ferrel cell etc)
    f = spint.interp1d(psi_coarse.lat, psi_coarse.sel(pfull=lev), axis=0, fill_value='extrapolate')
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
    get_edge_psi_nh('ap_2', sanity_check=True, thresh=120.)
    #get_edge_psi('ss_115.000', sanity_check=True, thresh=0., intdown=False)
    
    #get_edge_psi('ap_20', sanity_check=True, thresh=0.)
    #get_edge_psi('sn_2.000', sanity_check=True, thresh=0.)
    get_edge_psi_nh('full_qflux', lonin=[60.,150.], sanity_check=True, thresh=120.)



