'''8/08/2018 Correlate onset timing for (SCS, Indian, or BoB) monsoon onset with other ERA fields.
   Initially, try with 200 hPa relative vorticity'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
import pandas as pd


rcParams['figure.figsize'] = 15, 10
rcParams['font.size'] = 14

def correlate_onset(monsoon, var, lev=200, years=range(1979,2017)):
    
    plot_dir = '/scratch/rg419/plots/onset_variability/' + monsoon + '/' + var + '/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Choose monsoon to look at
    if monsoon == 'SCS':
        onsets = np.load('/scratch/rg419/python_scripts/python_bin_updates/physics_updates/climatology/isca_era_onsets_scsm.npy')
        onsets = np.array(onsets[1]) * 5. - 2.5
    elif monsoon == 'BOB':
        onsets = np.load('/scratch/rg419/python_scripts/python_bin_updates/physics_updates/climatology/era_onsets_bob_days.npy')
        onsets = np.array(onsets)
    elif monsoon == 'ISM':
        onsets = np.load('/scratch/rg419/python_scripts/python_bin_updates/physics_updates/climatology/era_onsets_ism_days.npy')
        onsets = np.array(onsets)
    else:
        print('Invalid input')
        return
    
    # Get mean onset date and standard deviation of onset date
    onset_stats = np.mean(np.array(onsets)), np.std(np.array(onsets))
    
    # Calculate difference of individual onsets from the mean
    onsets_diff = onsets - onset_stats[0]
    onsets_diff = xr.DataArray(onsets_diff, coords={'year': ('year', years)}, dims=['year'])
    
    # Load in ERA data, trim to only use 1979-2016, and resample in months
    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/sep_levs_' + var + '/era_' + var + '_' + str(lev) + '.nc')
    try:
        data = data[var].load().squeeze('level').loc['1979-01':'2016-12']  
    except:
        data = data[var].load().loc['1979-01':'2016-12']  
    
    data = data.resample(time='M').mean('time')
    
    # Use groupby to create a climatology
    data_clim = data.groupby('time.month').mean('time')
    
    # Calculate anomalies relative to the climatology for each month
    data_anom = data.groupby('time.month') - data_clim
    
    # Standard deviation is then rootmeansquare of anomalies
    data_std = np.sqrt((data_anom**2.).groupby('time.month').mean('time'))
    
    # Multiply ERA anomalies for each year by onset anomaly for that year, and average to get covariance
    covar = (onsets_diff * data_anom.groupby('time.year')).groupby('time.month').mean('time')
    
    # Correlation is the mean of the above product, normalised by the standard deviations of the onset dates and data
    correl = covar/onset_stats[1]/data_std
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    title = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    i=0
    for ax in axes:
        f1 = correl[i,:,:].plot.contourf(ax=ax, x='longitude',y='latitude', levels=np.arange(-1.,1.1,0.25), add_labels=False, add_colorbar=False)
        #correl[i,:,:].plot.contour(ax=ax, x='longitude',y='latitude', levels=np.arange(-100.4,100.,100.), colors='k', lw=2, add_labels=False)
        #correl[i,:,:].plot.contour(ax=ax, x='longitude',y='latitude', levels=np.arange(-99.6,100.,100.), colors='k', lw=2, add_labels=False)
        land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
        land = xr.open_dataset(land_mask)
        land.lsm[0,:,:].plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        ax.set_xticks(np.arange(0.,361.,90.))
        ax.set_yticks(np.arange(-90.,91.,30.))
        ax.set_title(title[i])
        ax.grid(True,linestyle=':')    
        i=i+1
    
    for ax in [ax1,ax5,ax9]:
        ax.set_ylabel('Latitude')
    
    for ax in [ax9,ax10,ax11,ax12]:
        ax.set_xlabel('Longitude')
    
    plt.subplots_adjust(left=0.06, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.2)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.15, aspect=60, shrink=0.75)
    cb1.set_label('Correlation coefficient')
    
    plt.savefig(plot_dir + 'corr_' + monsoon + '_' + var + '_' + str(lev) + '.pdf', format='pdf')
    
    plt.close()
    
    data.close()
    
    return correl
    

#for lev in range(50,301,50):
#    correlate_onset('BOB', 'vo', lev=lev)
#    correlate_onset('SCS', 'vo', lev=lev)
#    correlate_onset('ISM', 'vo', lev=lev)

#for lev in np.arange(50,1001,50):
    #print(lev)
    #correlate_onset('BOB', 'u', lev=lev)
    #correlate_onset('SCS', 'u', lev=lev)
    #correlate_onset('ISM', 'u', lev=lev)
    #correlate_onset('BOB', 'v', lev=lev)
    #correlate_onset('SCS', 'v', lev=lev)
    #correlate_onset('ISM', 'v', lev=lev)

correlate_onset('BOB', 'd', lev=150)
correlate_onset('SCS', 'd', lev=150)
correlate_onset('ISM', 'd', lev=150)
    