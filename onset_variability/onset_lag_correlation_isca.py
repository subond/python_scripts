'''9/08/2018 Correlate onset timing for some definition of Isca monsoon onset with other fields.'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
import pandas as pd
from climatology import isca_onsets
from data_handling_updates import isca_load_and_reshape

rcParams['figure.figsize'] = 15, 10
rcParams['font.size'] = 14

def correlate_onset_isca(run, var, onset_def='SCS', lonin=[110.,120.], land_mask=None, lev=150.):
    
    plot_dir = '/scratch/rg419/plots/onset_variability/isca_' + run + '/' + onset_def + '/' + var + '/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    
    # Load and reshape Isca data
    data = isca_load_and_reshape(run, var, months=[121,481])
    data.coords['month'] = (data.xofyear - 1) //6 + 1 
    data = data.groupby('month').mean(('xofyear'))
    
    # Create a climatology
    data_clim = data.mean('year_no')
    # Calculate anomalies relative to the climatology for each month
    data_anom = data - data_clim
    # Standard deviation is then rootmeansquare of anomalies
    data_std = np.sqrt((data_anom**2.).mean('year_no'))
    
    
    # Check if onset date dataset already exists
    try:
        savedir = '/scratch/rg419/onset_dates/'
        onsets = np.load(savedir + 'isca_onsets_' + run + '_' + onset_def + '_' + str(lonin[0]) + '_' + str(lonin[1]) + '.npy')
    # Otherwise produce it
    except:
        isca_onset(run, onset_def=onset_def, lonin=lonin)
        onsets = np.load(savedir + 'isca_onsets_' + run + '_' + onset_def + '_' + str(lonin[0]) + '_' + str(lonin[1]) + '.npy')
    onsets = np.array(onsets)
    # Get mean onset date and standard deviation of onset date
    onset_stats = np.mean(np.array(onsets)), np.std(np.array(onsets))
    # Calculate difference of individual onsets from the mean
    onsets_diff = onsets - onset_stats[0]
    onsets_diff = xr.DataArray(onsets_diff, coords={'year_no': ('year_no', data.year_no)}, dims=['year_no'])
    
    
    # Multiply ERA anomalies for each year by onset anomaly for that year, and average to get covariance
    covar = (onsets_diff * data_anom).mean('year_no')
    # Correlation is the mean of the above product, normalised by the standard deviations of the onset dates and data
    correl = covar/onset_stats[1]/data_std
    
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    title = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    i=0
    for ax in axes:
        f1 = correl.sel(pfull=lev)[i,:,:].plot.contourf(ax=ax, x='lon',y='lat', levels=np.arange(-1.,1.1,0.25), add_labels=False, add_colorbar=False)
        #correl[i,:,:].plot.contour(ax=ax, x='longitude',y='latitude', levels=np.arange(-100.4,100.,100.), colors='k', lw=2, add_labels=False)
        #correl[i,:,:].plot.contour(ax=ax, x='longitude',y='latitude', levels=np.arange(-99.6,100.,100.), colors='k', lw=2, add_labels=False)
        if not land_mask==None:
            land = xr.open_dataset(land_mask)
            land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
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
    
    plt.savefig(plot_dir + 'corr_' + onset_def + '_' + var + '_' + str(lev) + '_' + str(int(lonin[0])) + '_' + str(int(lonin[1])) + '.pdf', format='pdf')
    
    plt.close()
    
    data.close()
    
    return correl
    

correlate_onset_isca('half_shallow', 'ucomp', onset_def='SCS', lonin=[100.,140.], land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
correlate_onset_isca('half_shallow', 'ucomp', onset_def='SCS', lonin=[150.,160.], land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
correlate_onset_isca('half_shallow', 'ucomp', onset_def='SCS', lonin=[170.,180.], land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')

correlate_onset_isca('half_shallow', 'ucomp', onset_def='BOB', lonin=[100.,140.], land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
correlate_onset_isca('half_shallow', 'ucomp', onset_def='BOB', lonin=[150.,160.], land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
correlate_onset_isca('half_shallow', 'ucomp', onset_def='BOB', lonin=[170.,180.], land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')

