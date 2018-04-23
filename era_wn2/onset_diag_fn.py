# 4/01/2018 copied this from python_bin/physics/old plus added in the plotting function from clean_diags/onset_plot.py
#define monsoon onset as the pentad when rainfall of more than 8mm/day is observed
#read in precip data for a given run and locate onset pentad for each year at each gridpoint and produce pcolor plot

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams



name_temp = '/scratch/rg419/obs_and_reanalysis/datafiles/gpcp_1dd_v1.2_p1d.%04d%02d.nc'
names = [name_temp % (m,n) for m in range( 1997, 2015) for n in range(1,13) ]

#read data into xarray 
data = xr.open_mfdataset( names, chunks={'time': 30})
land_mask='/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'



def onset_diag_fn(rain):
    
    nlat = len(data.precip.lat)
    nlon = len(data.precip.lon)
    
    if len(rain.time)==366:
        pentad = np.repeat(np.arange(1., 74.), 5)
        pentad = np.insert(pentad, 10, 2)    
    else:
        pentad = np.repeat(np.arange(1., 74.), 5)
        
    rain = rain.assign_coords(pentad = ('time', pentad))
        
        #print year
        #print rain.sel(time=str(year)+'-12-31')
        #print rain.time.labels
        #print rain.get_index('time')
        #rain = rain.assign_coords( day = ('time', np.arange(1.,367.)))
        #rain = rain.drop(366, dim='time')
    
#    rain = rain.assign_coords( pentad = ('time', np.repeat( np.arange(1., 74.), 5) )  )
    
    rain = rain.groupby('pentad').mean(('time'))
    
    
    #mask values of precip less than 6mm/day
    tot_p_clim_masked = np.ma.masked_less(rain.values,8.)

    #locate first unmasked value along pentad axis
    onset_index = np.ma.notmasked_edges(tot_p_clim_masked, axis=0)
    
    #onset_index[0][0][:] = onset_index[0][0][:]//5 +1
    #onset_index[1][0][:] = onset_index[0][0][:]//5 +1
    
    onset = np.zeros((nlat, nlon))
    onset[:] = np.nan
    
    for i in range(0,len(onset_index[0][1])):
        onset[ onset_index[0][1][i], onset_index[0][2][i] ] = onset_index[0][0][i] + 1
        
    onset_pentad = xr.DataArray(onset, [('lat', rain.lat), ('lon', rain.lon)])
    
    return onset_pentad, onset_index#, onset_end




def onset_plot(year, land_file):
    
    rcParams['figure.figsize'] = 15, 7.5
    
    plot_dir = '/scratch/rg419/plots/era_wn2/onset_dates/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    onset_pentad, onset_index = onset_diag_fn(data.precip.sel(time=str(year)))

    land = xr.open_dataset( land_file)
   # land_plot = xr.DataArray(land.lsm.values, [('lat', onset_pentad.lat), ('lon', onset_pentad.lon)])
    
    #Start plotting
    ax = onset_pentad.plot.pcolormesh(x='lon', y='lat',levels=np.arange(18.,49.,1.), extend='min', add_colorbar=False, cmap='Greys')
    #for i in range(0,len(onset_index[0][1])):
        #if (onset_pentad.lon[onset_index[0][2][i]] > 40 and onset_pentad.lon[onset_index[0][2][i]] < 150 and 
        #    onset_pentad.lat[onset_index[0][1][i]] > 0 and onset_pentad.lat[onset_index[0][1][i]] < 50):
        #    if (onset_index[0][0][i]>18 and onset_index[0][0][i]<=39):
        #        plt.text(onset_pentad.lon[onset_index[0][2][i]],onset_pentad.lat[onset_index[0][1][i]],onset_index[0][0][i]+1,size=10,rotation='vertical')
        #    elif (onset_index[0][0][i]>=40 and onset_index[0][0][i]<47):
        #        plt.text(onset_pentad.lon[onset_index[0][2][i]],onset_pentad.lat[onset_index[0][1][i]],onset_index[0][0][i]+1,size=10,rotation='vertical')
    
    land.lsm[0,:,:].plot.contour(x='longitude', y='latitude',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    cb1=plt.colorbar(ax)
    cb1.set_label('Onset pentad')
    plt.ylim(0,50)
    plt.xlim(40,150)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.tight_layout()  
    plt.savefig( plot_dir + 'onset_pentad_' + str(year) + '.png')
    plt.clf()


for year in range(1997, 2015):
    print year
    onset_plot(year, land_mask)
    