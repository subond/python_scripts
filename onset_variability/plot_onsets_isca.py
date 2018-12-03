'''2/08/2018 Plot onset timing for SCS, Indian, and BoB monsoon onset, for comparison with each other and for validation of code against papers
(ERA-Interim u-850 based metrics)'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm

plot_dir = '/scratch/rg419/plots/onset_variability/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

def plot_onsets_isca(run, color_in='C0', lonin=[0.,360.]):
    
    savedir = '/scratch/rg419/onset_dates/'
    # Load in onset dates
    onsets_scsm = np.load(savedir + 'isca_onsets_' + run + '_SCS_' + str(lonin[0]) + '_' + str(lonin[1]) + '.npy')    
    #onsets_bob = np.load('/scratch/rg419/onset_dates/isca_onsets_' + run + '_BOB_0.0_360.0.npy')
    #onsets_ism = np.load('/scratch/rg419/onset_dates/isca_onsets_' + run + '_ISM_0.0_360.0.npy')
    
    # Take mean and standard deviation
    scsm_stats = np.mean(np.array(onsets_scsm)), np.std(np.array(onsets_scsm))
    
    years = range(30)
    # Decide whether monsoon is early, normal, or late
    def get_timing(onsets, years):
        onset_stats = np.mean(np.array(onsets)), np.std(np.array(onsets))
        print(onset_stats[1])
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
    #print(onsets_scsm.swap_dims({'year': 'timing'}).sel(timing='Normal').year)
    
    onsets_scsm = get_timing(onsets_scsm, years)
    #onsets_bob = get_timing(onsets_bob, years)
    #onsets_ism = get_timing(onsets_ism, years)
    
    onsets_scsm.plot(color=color_in, marker='x')
    #onsets_bob.plot(color=color_in, marker='x')
    #onsets_ism.plot(color='C2', marker='x')
    
    #plt.figure(2)
    #plt.plot(onsets_scsm, onsets_bob, 'x')
    
    #plt.figure(3)
    #plt.plot(onsets_scsm, onsets_ism, 'x')
    
    #plt.figure(4)
    #plt.plot(onsets_ism, onsets_bob, 'x')
    
    #plt.savefig(plot_dir + 'onset_timings_era.pdf', format='pdf')
    #plt.close()

plot_onsets_isca('mld_2.5')
plot_onsets_isca('mld_5', color_in='C1')
plot_onsets_isca('sn_1.000', color_in='C2')
plot_onsets_isca('mld_15', color_in='C3')
plot_onsets_isca('half_shallow', color_in='C4', lonin=[170.,180.])
plot_onsets_isca('half_shallow', color_in='C5', lonin=[150.,160.])
plot_onsets_isca('half_shallow', color_in='C6', lonin=[100.,140.])

plt.ylim(30,55)
plt.show()
