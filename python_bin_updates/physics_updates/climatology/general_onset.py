""" 3/08/2018 Use Wang and Linho 2002 relative pentad mean rainfall rate to calculate onset dates for isca data
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

plot_dir = '/scratch/rg419/plots/zonal_asym_runs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def get_onset(rain_rel): 
    rain_rel_masked = np.ma.masked_less(rain_rel.values, 5)  # Mask relative rainfall where it is less than 5mm/day
    onset_index = np.ma.notmasked_edges(rain_rel_masked, axis=0) # Look along pentad axis to find edges of mask
    onset = np.zeros((len(rain_rel.lat),len(rain_rel.lon))) # Create array of nans to load onsets into
    onset[:] = np.nan
    
    for i in range(0,len(onset_index[0][1])): # Extract onsets from mask edges
        onset[ onset_index[0][1][i], onset_index[0][2][i] ] = onset_index[0][0][i]+1
    onset_pentad = xr.DataArray(onset, [('lat', rain_rel.lat), ('lon', rain_rel.lon)]) # Make onset pentad into a dataarray
    return onset_pentad
        
         
def relative_rain(run, real_conts=False, land_mask_name=None):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    jan_precip = data.precipitation.sel(xofyear=range(1,7)).mean('xofyear')
    relative_rain = (data.precipitation - jan_precip)*86400.
    
    onset_pentad = get_onset(relative_rain)
    
    # Create figure with 2 subplots and plot up results
    rcParams['figure.figsize'] = 10, 3
    rcParams['font.size'] = 14
    fig, ax1 = plt.subplots(1, 1)
        
    f1 = onset_pentad.plot.pcolormesh(ax=ax1, x='lon',y='lat',  levels=np.arange(25.,56.), cmap='plasma_r', add_colorbar=False, add_labels=False)
    
    plt.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.18, hspace=0.1, wspace=0.2)
    
    cbar1 = plt.colorbar(f1, ax=ax1, ticks=np.arange(25.,56.,5.), fraction=0.05)    
    
    if real_conts:
        land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc'
    else:
        if land_mask_name==None:
            land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/' + run + '.nc'
        else:
            land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/' + land_mask_name + '.nc'
    
    land = xr.open_dataset(land_mask)
    land.land_mask.plot.contour(ax=ax1, x='lon', y='lat', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    ax1.grid(True,linestyle=':')
    ax1.set_ylim([0,50])
    ax1.set_xticks(np.arange(0.,361.,60.))
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    
    plt.savefig(plot_dir + 'rel_rain_onset_' + run + '.pdf', format='pdf')
    plt.close()
    



def relative_rain_cmap():
    '''10/10/2018 Load up CMAP pentad and monthly means, and use to evaluate Northern and Southern hemisphere monsoon onset using the Wang and LinHo criteria'''
    
    # Load up data
    data = xr.open_dataset('/disca/share/rg419/CMAP_precip.pentad.mean.nc')
    data_mon = xr.open_dataset('/disca/share/rg419/CMAP_precip.mon.mean.nc')
    
    # Add pentad as a coordinate and make a climatology
    data.coords['pentad'] = (('time'), np.tile(np.arange(1,74),38))
    precip_clim_mon = data_mon.precip.groupby('time.month').mean('time')
    precip_clim_pentad = data.precip.groupby('pentad').mean('time')
    
    # Calculate rainfal relative to January and July. For July re-order pentads so onset can be identified for Southern hemisphere monsoons too.
    rain_rel_jan = xr.DataArray((precip_clim_pentad - precip_clim_mon.sel(month=1)).values, [('pentad', np.arange(1,74)), ('lat', data.lat), ('lon', data.lon)])
    rain_rel_july = xr.DataArray((precip_clim_pentad - precip_clim_mon.sel(month=7)).values, [('pentad', (np.arange(1,74)-37)%73-36), ('lat', data.lat), ('lon', data.lon)])
    rain_rel_july = rain_rel_july.sortby('pentad')        
    
    onset_nh = get_onset(rain_rel_jan)
    onset_sh = get_onset(rain_rel_july) - 37.   # Onset is output as an index along the time axis. For July, where this has been rotated, therefore subtract 37 to correct
    
    # Create figure with 2 subplots and plot up results
    rcParams['figure.figsize'] = 10, 6
    rcParams['font.size'] = 14
    fig, (ax1,ax2) = plt.subplots(2, 1, sharex='col')
        
    f1 = onset_nh.plot.pcolormesh(ax=ax1, x='lon',y='lat',  levels=np.arange(20.,51.), cmap='plasma_r', add_colorbar=False, add_labels=False)
    f2 = onset_sh.plot.pcolormesh(ax=ax2, x='lon',y='lat', levels=np.arange(-15.,16.), cmap='plasma_r', add_colorbar=False, add_labels=False)
    
    plt.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.1, hspace=0.1, wspace=0.2)
    
    cbar1 = plt.colorbar(f1, ax=ax1, ticks=np.arange(20.,51.,5.), fraction=0.05)    
    cbar2 = plt.colorbar(f2, ax=ax2, ticks=np.arange(-15.,16.,5.), fraction=0.05)
    cbar2.ax.set_yticklabels((np.arange(-15,16,5)+73)%73)  # Label colorbar with pentad numbers corresponding to the onset pentad
    
    land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
    land = xr.open_dataset(land_mask)
    land.lsm[0,:,:].plot.contour(ax=ax1, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    land.lsm[0,:,:].plot.contour(ax=ax2, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    
    for ax in [ax1,ax2]:
        ax.grid(True,linestyle=':')
        ax.set_ylabel('Latitude')
        
    ax1.set_ylim([0,50])
    ax2.set_ylim([-50,0])
    ax2.set_xticks(np.arange(0.,361.,60.))
    ax2.set_xlabel('Longitude')
    
    plt.savefig(plot_dir + 'rel_rain_onset_cmap.pdf', format='pdf')
    plt.close()




def relative_rain_gpcp(nh=True, levels=np.arange(20.,51.,1)):
    # NOTE detrended GPCP is not really what you want - maybe fix this later but for now it's really just to keep the code work
    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/gpcp_detrendedpentads.nc')
    jan_precip = data.precip_clim.sel(xofyear=range(1,7)).mean('xofyear')
    july_precip = data.precip_clim.sel(xofyear=range(37,43)).mean('xofyear')
    
    if nh:
        relative_rain = (data.precip_clim - jan_precip)
    else:
        relative_rain = (data.precip_clim - july_precip)
        #print(relative_rain[36,:,:])
        relative_rain = xr.DataArray(relative_rain.values, [('xofyear', (data.xofyear-37)%73-36), ('lat', data.lat), ('lon', data.lon)])
        relative_rain = relative_rain.sortby('xofyear')
        #print(relative_rain[0,:,:])
        
    
    relative_rain_masked = np.ma.masked_less(relative_rain.values,5)
    
    #locate first unmasked value along pentad axis
    onset_index = np.ma.notmasked_edges(relative_rain_masked, axis=0)

    onset = np.zeros((180,360))
    onset[:] = np.nan
    #onset_end = np.zeros((64,128))
    #onset_end[:] = np.nan
    
    for i in range(0,len(onset_index[0][1])):
        if nh:
            onset[ onset_index[0][1][i], onset_index[0][2][i] ] = onset_index[0][0][i]+1
        else:
            onset[ onset_index[0][1][i], onset_index[0][2][i] ] = onset_index[0][0][i]-36
        #onset_end[ onset_index[1][1][i], onset_index[1][2][i] ] = onset_index[1][0][i]+1

    onset_pentad = xr.DataArray(onset, [('lat', data.lat), ('lon', data.lon % 360.)])
    onset_pentad = onset_pentad.sortby('lon')

    f = onset_pentad.plot.contourf(x='lon',y='lat', levels=levels, cmap='plasma_r', add_colorbar=False)
    
    cbar = plt.colorbar(f, ticks=levels[0::5])
    
    if not nh:
        cbar.ax.set_yticklabels((levels[0::5]+73)%73)
    
    land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
    land = xr.open_dataset(land_mask)
    land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
    #plt.xlim([40,180])
    #plt.xlim([0,100])
    if nh:
        plt.ylim([0,50])
    else:
        plt.ylim([-50,0])
    #plt.ylim([-50,0])
    
    rect = Rectangle((110, 5), 10, 10, edgecolor='m', facecolor='None', lw=2)
    ax = plt.gca()
    ax.add_patch(rect)

    rect = Rectangle((90, 5), 10, 10, edgecolor='m', facecolor='None', lw=2)
    ax = plt.gca()
    ax.add_patch(rect)
    
    rect = Rectangle((40, 5), 40, 10, edgecolor='m', facecolor='None', lw=2)
    ax = plt.gca()
    ax.add_patch(rect)
    
    if nh:
        plt.savefig(plot_dir + 'rel_rain_onset_gpcp.pdf', format='pdf')
    else:
        plt.savefig(plot_dir + 'rel_rain_onset_gpcp_sh.pdf', format='pdf')
    plt.close()




def relative_rain_gpcp_years(years=[1997]):
    
    plot_dir = '/scratch/rg419/plots/onset_variability/gpcp_years/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    # Load in all GPCP data
    name_temp = '/scratch/rg419/obs_and_reanalysis/datafiles/gpcp_1dd_v1.2_p1d.%04d%02d.nc'
    names = [name_temp % (m,n) for m in range( 1997, 2015) for n in range(1,13) ]
    data = xr.open_mfdataset( names, chunks={'time': 30})
    
    # Take mean over all jans
    #jan_precip = data.precip.sel(time=data['time.month']==1).mean('time')
    
    # Select requested year(s) and convert time to pentads
    for year in years:
        # Take mean over jan of the year wanted
        jan_precip = data.precip.sel(time=str(year)+'-01').mean('time')
        print(year)
        data_year = data.precip.sel(time=str(year))
        if len(data_year.time)==366:
            pentad = np.repeat(np.arange(1., 74.), 5)
            pentad = np.insert(pentad, 10, 2)    
        else:
            pentad = np.repeat(np.arange(1., 74.), 5)
        data_year = data_year.assign_coords(pentad = ('time', pentad))
        data_year = data_year.groupby('pentad').mean(('time'))
        
        relative_rain = (data_year - jan_precip)
        relative_rain_masked = np.ma.masked_less(relative_rain.values,5)
        #locate first unmasked value along pentad axis
        onset_index = np.ma.notmasked_edges(relative_rain_masked, axis=0)
        
        onset = np.zeros((180,360))
        onset[:] = np.nan
    
        for i in range(0,len(onset_index[0][1])):
            onset[ onset_index[0][1][i], onset_index[0][2][i] ] = onset_index[0][0][i]+1

        onset_pentad = xr.DataArray(onset, [('lat', data.lat), ('lon', data.lon)])
        
        onset_pentad.plot.contourf(x='lon',y='lat', levels=np.arange(15,45,2), extend='neither', cmap='Blues')
        land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
        land = xr.open_dataset(land_mask)
        land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        plt.xlim([40,180])
        #plt.xlim([0,100])
        plt.ylim([0,50])
        #plt.ylim([-50,0])
        
        rect = Rectangle((110, 5), 10, 10, edgecolor='m', facecolor='None', lw=2)
        ax = plt.gca()
        ax.add_patch(rect)

        rect = Rectangle((90, 5), 10, 10, edgecolor='m', facecolor='None', lw=2)
        ax = plt.gca()
        ax.add_patch(rect)
        
        rect = Rectangle((40, 5), 40, 10, edgecolor='m', facecolor='None', lw=2)
        ax = plt.gca()
        ax.add_patch(rect)
        
        plt.savefig(plot_dir + 'rel_rain_onset_gpcp_' + str(year) + '.pdf', format='pdf')
        plt.close()
        
    return


if __name__ == "__main__":
    
    #relative_rain('half_shallow')
    relative_rain('half_dry', land_mask_name='half_shallow')
    relative_rain('half_bright', land_mask_name='half_shallow')
    #relative_rain('half_shallow_5', land_mask_name='half_shallow')
    #relative_rain('half_shallow_10', land_mask_name='half_shallow')
    #relative_rain('half_nh_shallow')
    #relative_rain('q_shallow')
    #relative_rain('3q_shallow')
    
    #relative_rain_cmap()
    
    #relative_rain_gpcp_years(years=range(1997,2015))
    