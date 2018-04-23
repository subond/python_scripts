#define monsoon onset as the pentad when rainfall of more than 6mm/day is observed
#read in precip data for a given run and locate onset pentad for each year at each gridpoint and produce pcolor plot

import sys
sys.path.append('/scratch/rg419/workdir_moist/python_bin')
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from time_av_xr import load_year_xr

def onset_diag_fn(inp_fol,nproc):

    #set run name
    run_fol = inp_fol + '/np' + nproc
    year = 1
    rundata = load_year_xr(run_fol, [year])
    rundata.coords['pentad'] = (rundata.time // 5) + 1
    conv_p = rundata.convection_rain.groupby('pentad').mean(('time'))

    tot_p = xr.DataArray(np.zeros((72,64,128,5)), [('pentad', conv_p.pentad ), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    onset_pentad_yearly = xr.DataArray(np.zeros((64,128,5)), [('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])

    for year in range(2,7):
        rundata = load_year_xr(run_fol, [year])
        rundata.coords['pentad'] = (rundata.time // 5) + 1
        conv_p = rundata.convection_rain.groupby('pentad').mean(('time'))
        cond_p = rundata.condensation_rain.groupby('pentad').mean(('time'))
        tot_p[:,:,:,year-2] = (conv_p + cond_p)*86400

    tot_p_clim = tot_p.mean(('year'))

    #mask values of precip less than 9mm/day
    tot_p_clim_masked = np.ma.masked_less(tot_p_clim.values,9)

    #locate first unmasked value along pentad axis
    onset_index = np.ma.notmasked_edges(tot_p_clim_masked, axis=0)
    onset = np.zeros((64,128))
    onset[:] = np.nan
    for i in range(0,len(onset_index[0][1])):
        onset[onset_index[0][1][i],onset_index[0][2][i]] = onset_index[0][0][i]+1

    onset_pentad = xr.DataArray(onset, [('lat', rundata.lat), ('lon', rundata.lon)])

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/flat_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', rundata.lat), ('lon', rundata.lon)])
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
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/onset_pentad.png')
    plt.clf()

    return onset_pentad

