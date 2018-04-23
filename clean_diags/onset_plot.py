# Produce plots of the zonal mean terms in the momentum budget as a function of time for the aquaplanet simulations, showing the transition and change in balance of terms

from data_handling import time_means, season_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from physics import onset_diag_fn

rcParams['figure.figsize'] = 20, 10


def onset_plot(run, months, filename='plev_pentad', period_fac=1.):
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #data = time_means(run, months, filename=filename, timeav='pentad', period_fac=period_fac)
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    data['rain'] = (('xofyear','lat','lon'), (data.convection_rain + data.condensation_rain)*86400.)

#    onset_pentad, onset_index, onset_end = onset_diag_fn( data.rain )
    onset_pentad, onset_index = onset_diag_fn( data.rain )

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', onset_pentad.lat), ('lon', onset_pentad.lon)])
    #plot
    ax = onset_pentad.plot.pcolormesh(x='lon', y='lat',levels=np.arange(18.,49.,3.), extend='min', add_colorbar=False)
    for i in range(0,len(onset_index[0][1])):
        if (onset_pentad.lon[onset_index[0][2][i]] > 40 and onset_pentad.lon[onset_index[0][2][i]] < 150 and 
            onset_pentad.lat[onset_index[0][1][i]] > 0 and onset_pentad.lat[onset_index[0][1][i]] < 50 and 
            onset_index[0][0][i]>18 and onset_index[0][0][i]<47):
            plt.text(onset_pentad.lon[onset_index[0][2][i]],onset_pentad.lat[onset_index[0][1][i]],onset_index[0][0][i]+1,size=8,rotation='vertical')

    land.land_mask.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    cb1=plt.colorbar(ax)
    cb1.set_label('Onset pentad')
    plt.ylim(0,50)
    plt.xlim(40,150)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.tight_layout()  
    plt.savefig( plot_dir + 'onset_pentad.png')
    plt.clf()
    
    
    #check when precip decreases below mask again
    #ax = onset_end.plot.pcolormesh(x='lon', y='lat',levels=np.arange(18.,49.,3.), extend='min', add_colorbar=False)
    #for i in range(0,len(onset_index[0][1])):
    #    if (onset_end.lon[onset_index[1][2][i]] > 40 and onset_end.lon[onset_index[1][2][i]] < 150 and 
    #        onset_end.lat[onset_index[1][1][i]] > 0 and onset_end.lat[onset_index[1][1][i]] < 50 and 
    #        onset_index[1][1][i]>18 and onset_index[1][1][i]<47):
    #        plt.text(onset_pentad.lon[onset_index[1][2][i]],onset_pentad.lat[onset_index[1][1][i]],onset_index[1][0][i]+1,size=8,rotation='vertical')

    #land.land_mask.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #cb1=plt.colorbar(ax)
    #cb1.set_label('Onset end check')
    #plt.ylim(0,50)
    #plt.xlim(40,150)
    #plt.xlabel('Longitude')
    #plt.ylabel('Latitude')
    #plt.tight_layout()  
    #plt.savefig( plot_dir + 'onset_end_check.png')
    #plt.clf()
    
    #if lonrange[0] == 'all':   
    #    onset_date = np.mean(onset_pentad[37,:])
    #else:
    #    lonmin =  np.min([k for k, j in enumerate(onset_pentad.lon) if j >= lonrange[0] ]) 
    #    lonmax =  np.max([k for k, j in enumerate(onset_pentad.lon) if j <= lonrange[1] ])
    #    onset_date = np.mean(onset_pentad[37,lonmin:lonmax])
    #print onset_date
        
#onset_plot('ap_2', [121,481])
#onset_plot('ap_20_qflux', [121,481])
onset_plot('am_qflux', [121,481])
onset_plot('flat_qflux', [121,481])
#onset_plot('full_qflux', [121,481])
