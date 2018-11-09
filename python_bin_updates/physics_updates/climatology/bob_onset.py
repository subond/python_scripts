""" 2/08/2018 Function to calculate the Mao and Wu 2006 Bay of Bengal monsoon onset diagnostic.
NB - I am using pentads here, but the paper in question uses days. Will check against other papers, may need to switch to days
Takes in u field, and averages this between 5 and 15N, and 90-100 E if a lon dimension is given. 
Time must be in pentads
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from data_handling_updates import isca_load_and_reshape


def scsm_onset(u, pentaddim='xofyear', pdim='pfull', latdim='lat', londim='lon', print_onset=False):
    
    # Check if the input pentad name is in the xarray dimensions, rename to pentad if so
    if pentaddim in u.dims:
        u = u.rename(new_name_or_name_dict={pentaddim: 'pentad'})
    else:    # Otherwise cause a crash 
        u['pentaddim']
    if latdim in u.dims: # Same with lat
        u = u.rename(new_name_or_name_dict={latdim: 'lat'})
    else:    
        u['latdim']
    if londim in u.dims: # Same with lon
        u = u.rename(new_name_or_name_dict={londim: 'lon'})
    else:    
        u['londim']
        
    # Select 850 hPa pressure level, if multiple pressure levels are present
    if pdim in u.dims:
        u = u.rename(new_name_or_name_dict={pdim: 'pfull'})
        u = u.sel( pfull = 850.)
    
    coslat = np.cos(u.lat * np.pi/180.)
    sinlat = np.sin(u.lat * np.pi/180.)
    
    u_area_weighted = u * coslat # Weight u by coslat, to do area weighting
    
    # Get specified lats and lons, select pentad range to look at
    lats = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= 5. and u.lat[i] < 15.]
    lons = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= 90. and u.lon[i] < 100.]
    pentads = [u.pentad[i].values for i in range(len(u.pentad)) if u.pentad[i] >=  0]
        
    # Calculate the area mean u
    u_mean = u_area_weighted.sel(lon=lons).sel(lat=lats).sum(('lat','lon')) / (coslat.sel(lat=lats).sum('lat') * len(lons))
        
    for i in range(len(pentads)):  # For each pentad
        if ((u_mean.sel(pentad=pentads[i]).values) > 0.  # Check if u_mean is greater than zero
             #and (u_mean.sel(pentad=pentads[i:i+4]).mean('pentad') > 1.) # and if the u_mean over that pentad and next 3 is greater than 1
             and (np.sum(u_mean.sel(pentad=pentads[i:i+2]).values > 0.) >= 2.)): # and if u_mean is greater than zero in next 2 pentads (10 days)
            onset_pentad = pentads[i]  # If all that is true, that's your onset pentad
            if print_onset:
                print('Onset pentad: ', onset_pentad) # Print it
            return u_mean, onset_pentad # Return it, and u
    if print_onset:
        print('No onset') # If you loop the whole way through and there's no onset, i.e. code still hasn't returned an onset date and finished, print that if wanted
    return u_mean, None # And return u and an empty value

    

if __name__ == "__main__":
    
    u = isca_load_and_reshape('control_qflux', 'ucomp', months=[121,481])
    
    onset_isca=[]
    for year_no in u.year_no.values:
        print(year_no, ' of ', u.year_no.max().values)
        u_mean, onset_pentad = scsm_onset(u.sel(year_no=year_no))
        onset_isca.append(onset_pentad)
        
    
    filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_850.nc'
    data = xr.open_dataset(filename, chunks={'latitude': 100, 'longitude': 100})
    data = data.resample('D', dim='time', how='mean')
        
    def pentad_means_of_year(data, year):
        data_year = data.sel(time=str(year))
        if len(data_year.time)==366:
            pentad = np.repeat(np.arange(1., 74.), 5)
            pentad = np.insert(pentad, 10, 2)    
        else:
            pentad = np.repeat(np.arange(1., 74.), 5)
        data_year = data_year.assign_coords(pentad = ('time', pentad))
        data_year = data_year.groupby('pentad').mean(('time'))
        return data_year
        
    onset_era=[]    
    for year in range(1979,2017):
        print(year)
        u_year = pentad_means_of_year(data.u, year)
        u_era = scsm_onset(u_year, pdim='level', pentaddim='pentad', latdim='latitude', londim='longitude')
        onset_era.append(u_era[1])
        
    np.save('isca_era_onsets_bob',np.array([onset_isca, onset_era]))
    
    
    #data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/control_qflux_0.500.nc')
    #u_p5 = scsm_onset(data.ucomp)
    
    #data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/control_qflux_0.750.nc')
    #u_p75 = scsm_onset(data.ucomp)
    
    #data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/control_qflux.nc')
    #u_1 = scsm_onset(data.ucomp)
    
    #data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/control_qflux_1.250.nc')
    #u_1p25 = scsm_onset(data.ucomp)
    
    #data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/control_qflux_1.500.nc')
    #u_1p5 = scsm_onset(data.ucomp)
    
    #u_1[0].plot()
    #u_p5[0].plot()
    #u_p75[0].plot()
    #u_1p25[0].plot()
    #u_1p5[0].plot()
    #plt.show()