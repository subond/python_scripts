# Produce a plot of onset pentad for each grid box for the monsoon region for the topography experiment

#define monsoon onset as the pentad when rainfall of more than 6mm/day is observed
#read in precip data for a given run and locate onset pentad for each year at each gridpoint and produce pcolor plot

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import load_year_xr
from physics import onset_diag_fn

def onset_plot(inp_fol,years,lonrange=['all']):

    onset_pentad, onset_index = onset_diag_fn(inp_fol, years)

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/flat_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', onset_pentad.lat), ('lon', onset_pentad.lon)])
    #plot
    onset_pentad.plot.pcolormesh(x='lon', y='lat',levels=np.arange(34.,48.,1.))
    for i in range(0,len(onset_index[0][1])):
        if (onset_pentad.lon[onset_index[0][2][i]] > 60 and onset_pentad.lon[onset_index[0][2][i]] < 180 and 
            onset_pentad.lat[onset_index[0][1][i]] > 0 and onset_pentad.lat[onset_index[0][1][i]] < 30 and 
            onset_index[0][0][i]>5):
            plt.text(onset_pentad.lon[onset_index[0][2][i]],onset_pentad.lat[onset_index[0][1][i]],onset_index[0][0][i]+1,size=8,rotation='vertical')

    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    plt.ylim(0,30)
    plt.xlim(60,180)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Onset pentad')
    plt.savefig('/scratch/rg419/plots/mom_budg_work/onset_pentad_' + inp_fol + '.png')
    plt.clf()
    
    if lonrange[0] == 'all':   
        onset_date = np.mean(onset_pentad[37,:])
    else:
        lonmin =  np.min([k for k, j in enumerate(onset_pentad.lon) if j >= lonrange[0] ]) 
        lonmax =  np.max([k for k, j in enumerate(onset_pentad.lon) if j <= lonrange[1] ])
        onset_date = np.mean(onset_pentad[37,lonmin:lonmax])
    print onset_date
        
    

onset_plot('topo_10m',range(21,41),[60.,120.])
onset_plot('flat_10m',range(21,41),[60.,120.])
onset_plot('aquaplanet_2m',range(21,41))
onset_plot('aquaplanet_10m',range(21,41))
