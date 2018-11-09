"""
Locate the division between the Hadley cells
Option to only consider the max of the SH Hadley cell in both time dependendent and steady state cases
Returns latitude of division, and magnitude of overturning, with sign convention so that SH cell is positive
3/11/2017
Updated to take data rather than run name as an input 6/11/2017
"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from hadley_cell import mass_streamfunction
import scipy.interpolate as spint

# Define latitudes to interpolate psi onto
lats = np.arange(-45, 45.1, 0.1)


def get_edge_psi(data, lonin=[-1.,361.], sanity_check=False, lev=500., thresh=0., do_month_av=False, period_fac=1., intdown=True, nh=False):
    """Front end so same function can be used for both climatological and steady state data
    Inputs: run = run name
           lonin = longitudes to average over
           sanity_check = plot up overturning and latitudes of max and division to check
           lev = level to look at
           thresh = threshold to use for boundary. If = None a 0.05 of the maximum at this time is used
           do_month_av = option to take a monthly average and return this, if wanted
           period_fac = ratio of year length to Earth year
           intdown = direction to integrate for mass streamfunction
           nh = Only calculate for northern hemisphere """
    
    #Load data
    #if run=='era':
    #    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_v_alllevs.nc')
        # Take pentad means
    #    data.coords['xofyear'] = data.day_of_yr //5 + 1.  
    #    data = data.groupby('xofyear').mean(('day_of_yr'))
    #else:
    #    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
            
    # Choose lons to average over
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    # If data has a time dimension, use time series fn. Otherwise use steady state fn
    if 'xofyear' in data.coords:
        edge_loc, psi_max, psi_max_loc = get_edge_psi_time(data, lons, sanity_check=sanity_check, 
                             lev=lev, thresh=thresh, intdown=intdown, nh=nh)
                             
        # Take monthly average if wanted, return outputs
        if do_month_av:    
            psi_max.coords['month'] = np.mod( psi_max.xofyear -1, 72.*period_fac) //(6*period_fac) +1    
            psi_max_month = psi_max.groupby('month').mean(('xofyear'))
            edge_loc.coords['month'] = np.mod( edge_loc.xofyear -1, 72.*period_fac) //(6*period_fac) +1    
            edge_loc_month = edge_loc.groupby('month').mean(('xofyear'))
        
            return edge_loc, psi_max, edge_loc_month, psi_max_month
        else:
            return edge_loc, psi_max, psi_max_loc
            
    else:
        edge_loc, psi_max, psi_max_loc = get_edge_psi_ss(data, lons, sanity_check=sanity_check, 
                             lev=lev, thresh=thresh, intdown=intdown, nh=nh)
        return edge_loc, psi_max, psi_max_loc
                             



def get_edge_psi_time(data, lons, sanity_check=False, lev=500., thresh=0., intdown=True, nh=False):
    """Get cell edge for climatological data    
       Inputs: data = input data file
           lons = longitudes to average over
           sanity_check = plot up overturning and latitudes of max and division to check
           lev = level to look at
           thresh = threshold to use for boundary
           intdown = direction to integrate for mass streamfunction
           nh = Only calculate for northern hemisphere """
    
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
        
    elif thresh <= 1.:
        # If a fraction of the max streamfunction is specified, use this to define the mask
        if nh:
            psi_max = psi.min('lat')
            thresh_frac = psi_max.values*-thresh
        else:
            psi_max = np.abs(psi).max('lat')
            thresh_frac = psi_max.values*thresh
        
        for i in range(psi.xofyear.shape[0]):
            psi_mask[np.abs(psi[:,i]).values <= thresh_frac[i], i] = np.nan
                        
    else:
        psi_mask[np.abs(psi).values <= thresh] = np.nan

    psi_mask = xr.DataArray( psi_mask, psi.coords)
    
    #plt.figure(1)
    #psi.plot.contourf(x='xofyear', y='lat')
    #plt.figure(2)
    #psi_mask.plot.contourf(x='xofyear', y='lat')
    #plt.show()
    
    # Use mask to mask psi at 500hPa 
    psi_masked = psi * psi_mask
    
    # Get a list of times where abs(psi) has values above the threshold
    times_defd = []
    for i in range(0, len(psi_masked.xofyear)):
        if np.any(np.isfinite(psi_masked[:,i])):
            times_defd.append(np.float(psi_masked.xofyear[i]))
    
    # Reduce psi to only include times where it goes above the threshold
    if thresh ==0.:
        psi_red = psi.sel(xofyear=times_defd)
    else:
        psi_red = psi_masked.sel(xofyear=times_defd)
    
    
    psi_mask_red = np.nan_to_num(psi_mask.sel(xofyear=times_defd))
    
    if nh:
        # If only looking at Northern Hemisphere, now take psi_red and look for min instead of taking abs and finding max    
        psi_max = psi_red.min('lat')
        psi_max_loc = psi_red.lat.values[psi_red.argmin('lat').values]
        psi_max_loc = xr.DataArray(psi_max_loc, coords=[psi_red.xofyear], dims=['xofyear'])
        ispos = []
        #Create a list of which times the min of psi_red is positive or negative
        for i in range(0, len(psi_max_loc)):
            ispos.append(psi_red.values[psi_red.argmin('lat').values[i],i] >= 0.)
        # Multiply by -1 to give a positive sign for the SH cell
        psi_max = -1.*psi_max
    else:
        psi_max = np.abs(psi_red).max('lat')
        psi_max_loc = psi_red.lat.values[np.abs(psi_red).argmax('lat').values]
        psi_max_loc = xr.DataArray(psi_max_loc, coords=[psi_red.xofyear], dims=['xofyear'])
        # Create list of which times the max of abs(psi_red) is positive vs negative
        ispos = []
        for i in range(0, len(psi_max_loc)):
            ispos.append(psi_red.values[np.abs(psi_red).argmax('lat').values[i],i] >= 0.)
        # Use this first to get the sign of psi_max right
        #NB this line not working properly on 26/01/2018. Use caution when this option is used.
        psi_max[ispos] = -1.*psi_max[ispos]
    
    
    # Now locate cell edge finally!
    edge_loc = np.zeros([len(psi_max_loc),])
    
    for i in range(0, len(psi_max_loc)):
        # Locate places where the psi mask changes sign    
        edges_are_at = psi_red.lat[np.where(psi_mask_red[:-1,i] != psi_mask_red[1:,i])[0]]
        try:
        # Use a try statement to avoid errors where edge is poorly defined
            if ispos[i]:
                # If the overturning is positive, look for the first edge to the south of the max
                # If only looking at the Northern hemisphere don't need this
                if nh:
                    edge_loc[i] = np.nan
                else:    
                    edge_loc[i] = np.max(edges_are_at[edges_are_at < psi_max_loc[i]])
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
        psi_masked.plot.contourf(levels=np.arange(-500.,510.,50.))
        psi_max_loc.plot()
        edge_loc.plot()
        
        plt.figure(2)
        plt.plot(edge_loc, psi_max, 'x')
        
        #plt.figure(3)
        #plt.plot(psi_max_loc.xofyear, psi_max_loc, 'x')
        plt.show()
    
    return edge_loc, psi_max, psi_max_loc






def get_edge_psi_ss(data, lons, sanity_check=False, lev=500., thresh=0., intdown=True, nh=False):
    """Get cell edge for steady state data
       Inputs: data = input data file
           lons = longitudes to average over
           sanity_check = plot up overturning and latitudes of max and division to check
           lev = level to look at
           thresh = threshold to use for boundary
           intdown = direction to integrate for mass streamfunction
           nh = Only calculate for northern hemisphere """
    
    # Calculate mass streamfunction
    psi_coarse = mass_streamfunction(data, a=6376.0e3, dp_in=50., lons=lons, intdown=intdown)
    psi_coarse /= 1.e9
    
    # Interpolate between 35S and 35N (Should help to avoid mistakes with Ferrel cell etc)
    f = spint.interp1d(psi_coarse.lat, psi_coarse.sel(pfull=lev), axis=0, fill_value='extrapolate')
    psi = f(lats)
    psi = xr.DataArray(psi, coords=[lats], dims=['lat'])
    
    # Create a mask that puts NaNs where abs(psi) is below a threshold
    psi_mask = np.ones(psi.values.shape)
    
    if thresh == 0.:
        psi_mask[psi.values <= thresh] = np.nan
    
    elif thresh == None:
        if nh:
            psi_max = psi.min('lat')
            thresh_frac = psi_max.values*-0.05
        else:
            psi_max = np.abs(psi).max('lat')
            thresh_frac = psi_max.values*0.05
        psi_mask[np.abs(psi).values <= thresh_frac] = np.nan
            
    else:
        psi_mask[np.abs(psi).values <= thresh] = np.nan
    
    psi_mask = xr.DataArray( psi_mask, psi.coords)
    
    # Use mask to mask psi at 500hPa 
    psi_masked = psi * psi_mask
    
    # Reduce psi to only where it goes above the threshold
    if thresh ==0.:
        psi_red = psi    
    else:
        psi_red = psi_masked
    
    psi_mask_red = np.nan_to_num(psi_mask)
    
    if nh:
        psi_max = psi_red.min('lat')
        psi_max_loc = psi_red.lat.values[psi_red.argmin('lat').values]
        # Create list of which times the max of abs(psi_rad) is positive vs negative
        ispos = psi_red.values[psi_red.argmin('lat').values] >= 0.
        
    else:
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

    return edge_loc, psi_max, psi_max_loc






#if __name__ == "__main__":
    #data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_100.000.nc')
    #vars = get_edge_psi(data, sanity_check=True, nh=True)
    #print vars
    #get_edge_psi('ap_2', sanity_check=True, thresh=None, nh=True)
    
    #get_edge_psi('ss_100.000', sanity_check=True, thresh=None, nh=False)
    #get_edge_psi('ss_100.000', sanity_check=True, thresh=0., nh=False)
    
    #get_edge_psi('ap_2', sanity_check=True, thresh=120.)
    #get_edge_psi('sn_2.000', sanity_check=True, thresh=0.)
    #get_edge_psi('full_qflux', lonin=[60.,150.], sanity_check=True, thresh=120.)



