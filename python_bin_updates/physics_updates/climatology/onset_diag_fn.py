# 4/01/2018 copied this from python_bin/physics/old plus added in the plotting function from clean_diags/onset_plot.py
#define monsoon onset as the pentad when rainfall of more than 8mm/day is observed
#read in precip data for a given run and locate onset pentad for each year at each gridpoint and produce pcolor plot

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams


def onset_diag_fn(rain):

    #mask values of precip less than 8mm/day
    tot_p_clim_masked = np.ma.masked_less(rain.values,6)

    #locate first unmasked value along pentad axis
    onset_index = np.ma.notmasked_edges(tot_p_clim_masked, axis=0)

    onset = np.zeros((64,128))
    onset[:] = np.nan
    #onset_end = np.zeros((64,128))
    #onset_end[:] = np.nan
    
    for i in range(0,len(onset_index[0][1])):
        onset[ onset_index[0][1][i], onset_index[0][2][i] ] = onset_index[0][0][i]+1
        #onset_end[ onset_index[1][1][i], onset_index[1][2][i] ] = onset_index[1][0][i]+1

    onset_pentad = xr.DataArray(onset, [('lat', rain.lat), ('lon', rain.lon)])
    #onset_end = xr.DataArray(onset_end, [('lat', rain.lat), ('lon', rain.lon)])

    return onset_pentad, onset_index#, onset_end




def onset_plot(run, land_file):
    
    rcParams['figure.figsize'] = 15, 7.5
    
    plot_dir = '/scratch/rg419/plots/climatology/onset_dates/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    # Open dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    try:
        precip = data.precipitation*86400.
    except:
        precip = (data.convection_rain + data.condensation_rain)*86400.
    
    onset_pentad, onset_index = onset_diag_fn(precip)

    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', onset_pentad.lat), ('lon', onset_pentad.lon)])
    
    #Start plotting
    ax = onset_pentad.plot.pcolormesh(x='lon', y='lat',levels=np.arange(18.,49.,1.), extend='min', add_colorbar=False, cmap='Greys')
    for i in range(0,len(onset_index[0][1])):
        if (onset_pentad.lon[onset_index[0][2][i]] > 40 and onset_pentad.lon[onset_index[0][2][i]] < 150 and 
            onset_pentad.lat[onset_index[0][1][i]] > 0 and onset_pentad.lat[onset_index[0][1][i]] < 50):
            if (onset_index[0][0][i]>18 and onset_index[0][0][i]<=39):
                plt.text(onset_pentad.lon[onset_index[0][2][i]],onset_pentad.lat[onset_index[0][1][i]],onset_index[0][0][i]+1,size=10,rotation='vertical')
            elif (onset_index[0][0][i]>=40 and onset_index[0][0][i]<47):
                plt.text(onset_pentad.lon[onset_index[0][2][i]],onset_pentad.lat[onset_index[0][1][i]],onset_index[0][0][i]+1,size=10,rotation='vertical', color='w')

    land.land_mask.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    cb1=plt.colorbar(ax)
    cb1.set_label('Onset pentad')
    plt.ylim(0,50)
    plt.xlim(40,150)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.tight_layout()  
    plt.savefig( plot_dir + 'onset_pentad_ ' + run+ '.png')
    plt.clf()




def onset_diff_plot(run_1, run_2, plot_title, land_file_1=None, land_file_2=None):
    
    rcParams['figure.figsize'] = 15, 7.5
    
    plot_dir = '/scratch/rg419/plots/climatology/onset_dates/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    # Open dataset
    data_1 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run_1+'.nc')
    data_2 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run_2+'.nc')
    
    try:
        data_1['precipitation'] = (data_1.precipitation*86400.)
    except:
        data_1['precipitation'] = ((data_1.convection_rain + data_1.condensation_rain)*86400.)
    
    try:
        data_2['precipitation'] = (data_2.precipitation*86400.)
    except:
        data_2['precipitation'] = ((data_2.convection_rain + data_2.condensation_rain)*86400.)
    
    onset_pentad_1, onset_index_1 = onset_diag_fn(data_1.precipitation)
    onset_pentad_2, onset_index_2 = onset_diag_fn(data_2.precipitation)
    
    onset_pentad_diff = onset_pentad_1 - onset_pentad_2
    
    land = xr.open_dataset( land_file_1)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', onset_pentad_1.lat), ('lon', onset_pentad_1.lon)])
    
    #onset_pentad_diff.plot.contourf()#levels=np.arange(18.,48.))
    #land_plot.plot.contour(colors='k')
    #plt.show()
    
    #Start plotting
    ax = onset_pentad_diff.plot.pcolormesh(x='lon', y='lat', levels=np.arange(-10.,10.5), extend='both', add_colorbar=False)
    #for i in range(0,len(onset_index[0][1])):
    #    if (onset_pentad.lon[onset_index[0][2][i]] > 40 and onset_pentad.lon[onset_index[0][2][i]] < 150 and 
    #        onset_pentad.lat[onset_index[0][1][i]] > 0 and onset_pentad.lat[onset_index[0][1][i]] < 50):
    #        if (onset_index[0][0][i]>18 and onset_index[0][0][i]<=39):
    #            plt.text(onset_pentad.lon[onset_index[0][2][i]],onset_pentad.lat[onset_index[0][1][i]],onset_index[0][0][i]+1,size=10,rotation='vertical')
    #        elif (onset_index[0][0][i]>=40 and onset_index[0][0][i]<47):
    #            plt.text(onset_pentad.lon[onset_index[0][2][i]],onset_pentad.lat[onset_index[0][1][i]],onset_index[0][0][i]+1,size=10,rotation='vertical', color='w')

    land.land_mask.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    cb1=plt.colorbar(ax)
    cb1.set_label('Onset pentad')
    plt.ylim(0,50)
    plt.xlim(40,150)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.tight_layout()  
    plt.savefig( plot_dir + 'onset_pentad_diff_' + plot_title + '.png')
    plt.clf()
    
    

if __name__ == "__main__":


    #onset_plot('idealised_1hill', '/scratch/rg419/Experiments/wavenumber_2/input/1hill.nc')
    #onset_plot('idealised_2hill', '/scratch/rg419/Experiments/wavenumber_2/input/2hill.nc')
    #onset_diff_plot('idealised_2hill', 'idealised_1hill', '2hill_1hill_diff',
    #         land_file_1='/scratch/rg419/Experiments/wavenumber_2/input/1hill.nc', land_file_2='/scratch/rg419/Experiments/wavenumber_2/input/2hill.nc')
    
    onset_plot('control_qflux', '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    onset_plot('no_americas', '/scratch/rg419/python_scripts/land_era/land_era_no_america.nc')
    onset_diff_plot('control_qflux', 'no_americas', 'full_no_am_diff',
             land_file_1='/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
             
    onset_plot('no_TIP', '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    onset_diff_plot('control_qflux', 'no_TIP', 'full_no_TIP_diff',
             land_file_1='/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    
    
    
    #onset_plot('full_qflux', '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')

    #onset_plot('idealised_land_sh_topo', '/scratch/rg419/Experiments/idealised_land/input/land_extra_topo.nc')
    #onset_plot('idealised_land_sh_topo_30m', '/scratch/rg419/Experiments/idealised_land/input/land_extra_topo.nc')

#onset_plot('idealised_land_ea_only', '/scratch/rg419/Experiments/idealised_land/input/land_EA.nc')
#onset_plot('idealised_land_in_only', '/scratch/rg419/Experiments/idealised_land/input/land_india.nc')
#onset_plot('idealised_land_no_af', '/scratch/rg419/Experiments/idealised_land/input/land_asia_tibet.nc')
#onset_plot('idealised_land_no_am', '/scratch/rg419/Experiments/idealised_land/input/land_aus_asia_tibet.nc')
#onset_plot('idealised_land_no_aus', '/scratch/rg419/Experiments/idealised_land/input/land_africa_asia_tibet.nc')
#onset_plot('idealised_land_no_sa', '/scratch/rg419/Experiments/idealised_land/input/land_EA_tibet.nc')
    