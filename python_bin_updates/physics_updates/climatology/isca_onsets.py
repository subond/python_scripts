""" 09/08/2018 Plot onset date for each year for a given lon range and definition
"""

import numpy as np
import xarray as xr
import sh
from data_handling_updates import isca_load_and_reshape


savedir = '/scratch/rg419/onset_dates/'
mkdir = sh.mkdir.bake('-p')
mkdir(savedir)
    
def scsm_onset(u, print_onset=False, lonin=[110.,120.]):
    
    u = u.sel( pfull = 850.)
        
    coslat = np.cos(u.lat * np.pi/180.)
    sinlat = np.sin(u.lat * np.pi/180.)
    
    u_area_weighted = u * coslat # Weight u by coslat, to do area weighting
    
    # Get specified lats and lons, select pentad range to look at
    lats = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= 5. and u.lat[i] <= 15.]
    lons = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= lonin[0] and u.lon[i] <= lonin[1]]
    pentads = [u.xofyear[i].values for i in range(len(u.xofyear)) if u.xofyear[i] >=  24]
    
    # Calculate the area mean u
    u_mean = u_area_weighted.sel(lon=lons).sel(lat=lats).sum(('lat','lon')) / (coslat.sel(lat=lats).sum('lat') * len(lons))
    
    for i in range(len(pentads)):  # For each pentad
        if ((u_mean.sel(xofyear=pentads[i]).values) > 0.  # Check if u_mean is greater than zero
             and (all(u_mean.sel(xofyear=pentads[i:i+4]).cumsum('xofyear').values/np.arange(1.,5.) >= 1.)) # and if the cumulative mean over that pentad and next 3 is greater than 1
             and (np.sum(u_mean.sel(xofyear=pentads[i:i+4]).values > 0.) >= 3.)): # and if u_mean is greater than zero in at least 3 out of 4 pentads
            onset_pentad = pentads[i]  # If all that is true, that's your onset pentad
            if print_onset:
                print('Onset pentad: ', onset_pentad) # Print it
            return u_mean, onset_pentad # Return it, and u
    if print_onset:
        print('No onset') # If you loop the whole way through and there's no onset, i.e. code still hasn't returned an onset date and finished, print that if wanted
    return u_mean, None # And return u and an empty value
    


def bob_onset(u, print_onset=False, lonin=[90.,100.]):
    
    u = u.sel( pfull = 850.)
        
    coslat = np.cos(u.lat * np.pi/180.)
    sinlat = np.sin(u.lat * np.pi/180.)
    
    u_area_weighted = u * coslat # Weight u by coslat, to do area weighting
    
    # Get specified lats and lons, select pentad range to look at
    lats = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= 5. and u.lat[i] <= 15.]
    lons = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= lonin[0] and u.lon[i] <= lonin[1]]
    pentads = u.xofyear.values
    
    # Calculate the area mean u
    u_mean = u_area_weighted.sel(lon=lons).sel(lat=lats).sum(('lat','lon')) / (coslat.sel(lat=lats).sum('lat') * len(lons))

    for i in range(len(pentads)):  # For each pentad
        if (all(u_mean.sel(xofyear=pentads[i:i+2]).values > 0)
            and (u_mean.sel(xofyear=pentads[i-1]).values < 0)): # Check if u_mean changes from negative to positive, and if it stays so for the next 10 days
            onset_pentad = pentads[i]  # If all that is true, that's your onset pentad
            if print_onset:
                print('Onset pentad: ', onset_pentad) # Print it
            return u_mean, onset_pentad # Return it, and u
    if print_onset:
        print('No onset') # If you loop the whole way through and there's no onset, i.e. code still hasn't returned an onset date and finished, print that if wanted
    return u_mean, None # And return u and an empty value



def ism_onset(u, print_onset=False, lonin=[40.,80.]):
    
    u = u.sel( pfull = 850.)
        
    coslat = np.cos(u.lat * np.pi/180.)
    sinlat = np.sin(u.lat * np.pi/180.)
    
    u_area_weighted = u * coslat # Weight u by coslat, to do area weighting
    
    # Get specified lats and lons, select pentad range to look at
    lats = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= 5. and u.lat[i] <= 15.]
    lons = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= lonin[0] and u.lon[i] <= lonin[1]]
    pentads = u.xofyear.values
    
    # Calculate the area mean u
    u_mean = u_area_weighted.sel(lon=lons).sel(lat=lats).sum(('lat','lon')) / (coslat.sel(lat=lats).sum('lat') * len(lons))
    
    for i in range(len(pentads)):  # For each pentad
        if (all(u_mean.sel(xofyear=pentads[i:i+2]).values > 6.2)):  # Check if u_mean is greater than 6.2m/s on that pentad and the next one
            onset_pentad = pentads[i]  # If that is true, that's your onset pentad
            if print_onset:
                print('Onset pentad: ', onset_pentad) # Print it
            return u_mean, onset_pentad # Return it, and u
    if print_onset:
        print('No onset') # If you loop the whole way through and there's no onset, i.e. code still hasn't returned an onset date and finished, print that if wanted
    return u_mean, None # And return u and an empty value
    
    

def isca_onset(run, lonin=[110.,120.], onset_def='SCS', print_onset=False):
    
    # Open all files and reshape
    u = isca_load_and_reshape(run, 'ucomp', months=[121,481])
    
    #Create empty list for onset dates
    onset_isca=[]
    
    # Based on 'def' call the relevant onset definition over the given region
    for year_no in u.year_no.values:
        print(year_no, ' of ', u.year_no.max().values)
        
        if onset_def=='SCS':
            u_mean, onset_pentad = scsm_onset(u.sel(year_no=year_no), print_onset=print_onset, lonin=lonin)
        elif onset_def=='BOB':
            u_mean, onset_pentad = bob_onset(u.sel(year_no=year_no), print_onset=print_onset, lonin=lonin)
        elif onset_def=='ISM':
            u_mean, onset_pentad = ism_onset(u.sel(year_no=year_no), print_onset=print_onset, lonin=lonin)
            
        onset_isca.append(onset_pentad)
    
    np.save(savedir + 'isca_onsets_' + run + '_' + onset_def + '_' + str(lonin[0]) + '_' + str(lonin[1]), np.array(onset_isca))
    return onset_isca


if __name__ == "__main__":
    
    #isca_onset('sn_1.000', lonin=[0., 360.])
    #isca_onset('sn_1.000', lonin=[0., 360.], onset_def='BOB')
    #isca_onset('sn_1.000', lonin=[0., 360.], onset_def='ISM')
    
    for run in [ 'mld_5', 'mld_15']:
        isca_onset(run, lonin=[0., 360.])
        isca_onset(run, lonin=[0., 360.], onset_def='BOB')
        isca_onset(run, lonin=[0., 360.], onset_def='ISM')
    
    
    #isca_onset('half_shallow', lonin=[170.,180.])
    #isca_onset('half_shallow', lonin=[150.,160.])
    #isca_onset('half_shallow', lonin=[100.,140.])
        
    #isca_onset('half_shallow', onset_def='BOB', lonin=[170.,180.])
    #isca_onset('half_shallow', onset_def='BOB', lonin=[150.,160.])
    #isca_onset('half_shallow', onset_def='BOB', lonin=[100.,140.])
    
    #isca_onset('half_shallow', onset_def='ISM', lonin=[170.,180.])
    #isca_onset('half_shallow', onset_def='ISM', lonin=[150.,160.])
    #isca_onset('half_shallow', onset_def='ISM', lonin=[100.,140.])
    
    
    