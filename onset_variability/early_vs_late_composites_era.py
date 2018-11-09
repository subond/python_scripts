'''25/09/2018 Identify early/late/normal onset timing for the different monsoons, and make a climatology of fields for these years using ERA data'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
import pandas as pd
import datetime

rcParams['figure.figsize'] = 15, 10
rcParams['font.size'] = 14

# Load in onset days
# Find years in which onset is early or late and make a list
# Select these years from the chosen dataset, average together and plot month by month data


def early_vs_late_onset(var, levs=[200], years=range(1979,2017), levels=np.arange(-6.,6.1,0.5)):
    # Load in onset dates
    onsets_scsm = np.load('/scratch/rg419/python_scripts/onset_variability/era_onsets_scsm.npy')
    onsets_bob  = np.load('/scratch/rg419/python_scripts/onset_variability/era_onsets_bob_days.npy')
    onsets_ism  = np.load('/scratch/rg419/python_scripts/onset_variability/era_onsets_ism_days.npy')
    onsets_scsm = onsets_scsm * 5. - 2.5
        
    # Take mean and standard deviation
    scsm_stats = np.mean(np.array(onsets_scsm)), np.std(np.array(onsets_scsm))
    bob_stats = np.mean(np.array(onsets_bob)), np.std(np.array(onsets_bob))
    ism_stats = np.mean(np.array(onsets_ism)), np.std(np.array(onsets_ism))
    
    # Decide whether monsoon is early, normal, or late
    def get_timing(onsets, years):
        onset_stats = np.mean(np.array(onsets)), np.std(np.array(onsets))
        #print(onset_stats[1])
        enl = []
        enl_no = []
        i=0
        for year in years:
            if onsets[i] > round(onset_stats[0] + 0.75*onset_stats[1]):
                enl.append('Late')
                enl_no.append(2)
            elif onsets[i] < round(onset_stats[0] - 0.75*onset_stats[1]):
                enl.append('Early')
                enl_no.append(0)
            else:
                enl.append('Normal')
                enl_no.append(1)
            
            i=i+1    
        onsets_out = xr.DataArray(onsets, coords={'timing': ('year', enl), 'timing_no': ('year', enl_no), 'year': ('year', years)}, dims=['year'])
        return onsets_out
        
    onsets_scsm = get_timing(onsets_scsm, years)
    onsets_bob = get_timing(onsets_bob, years)
    onsets_ism = get_timing(onsets_ism, years)
    
    # Load in ERA data, trim to only use 1979-2016, and resample in months
    name_temp = '/disca/share/reanalysis_links/monthly_era/' + var + '/interim_' + var + '_%04d.nc'
    names = [name_temp % m for m in years  ]
    
    data = xr.open_mfdataset(names)
    data = data[var].load().loc['1979-01':'2016-12']
    
    # Add a coordinate with the early/late/normal timing    
    data.coords['timing_scsm'] = (('time'), np.repeat(onsets_scsm.timing.values,12))
    data.coords['timing_bob'] = (('time'), np.repeat(onsets_bob.timing.values,12))
    data.coords['timing_ism'] = (('time'), np.repeat(onsets_ism.timing.values,12))    
    
    def plot_by_timing(timing_coord, timing, levels=levels):
        # Mask so only chosen timing remains, then groupby month and average
        data_timing = data.where(data[timing_coord]==timing, drop=True).groupby('time.month').mean('time') - data.groupby('time.month').mean('time')
        for lev in levs:
            try:
                data_timing = data_timing.sel(level=lev)
            except:
                data_timing = data_timing
                
            # Start figure with 12 subplots
            fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
            axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
            title = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
            i=0
            for ax in axes:
                f1 = data_timing[i,:,:].plot.contourf(ax=ax, x='longitude',y='latitude', add_labels=False, add_colorbar=False, extend='both',levels=levels)
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
        cb1.set_label(var)
        
        plot_dir = '/scratch/rg419/plots/onset_variability_new/early_vs_late/' + timing_coord + '/' + var + '/'
        mkdir = sh.mkdir.bake('-p')
        mkdir(plot_dir)
        
        plt.savefig(plot_dir + 'era_' + timing + '_' + var + '_' + str(lev) + '.pdf', format='pdf')
        plt.close()
    
    for monsoon in ['timing_scsm', 'timing_bob', 'timing_ism']:
        plot_by_timing(monsoon, 'Early', levels=levels)
        plot_by_timing(monsoon, 'Normal', levels=levels)
        plot_by_timing(monsoon, 'Late', levels=levels)
    
    
    data.close()
    

var_list = [['vo',np.arange(-2.e-5,2.1e-5,2.e-6)],
        ['u',np.arange(-5.,5.1,0.5)],
        ['v',np.arange(-5.,5.1,0.5)],
        ['t',np.arange(-1.2,1.21,0.2)],
        ['q',np.arange(-1.e-5,1.e-5,1.e-6)],
        ['w',np.arange(-0.03,0.031,0.005)],
        ['z',np.arange(-800.,801.,50.)]]

for var in var_list:
    early_vs_late_onset(var[0], levels=var[1], levs=[200])
    early_vs_late_onset(var[0], levels=var[1], levs=[850])