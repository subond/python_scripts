""" 14/06/2018 Plot histograms of the onset dates in isca and in era data and save
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from climatology import scsm_onset
from data_handling_updates import isca_load_and_reshape

rcParams['figure.figsize'] = 7, 5
rcParams['font.size'] = 14

plot_dir = '/scratch/rg419/plots/scs_monsoon/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def get_onsets(run, months=[121,481], period_fac=1.):
    
    year_starts = np.arange(months[0], months[1], 12)
    onset=[]
    for year_start in year_starts:
        name_temp = '/disca/share/rg419/Data_moist/' + run + '/run%04d/plev_pentad.nc'
        names = [name_temp % m for m in range(year_start, year_start + 12)]
        #read data into xarray 
        data = xr.open_mfdataset( names,
            decode_times=False,  # no calendar so tell netcdf lib
            # choose how data will be broken down into manageable chunks.
            chunks={'time': 30})
        data.coords['xofyear'] = np.mod( data.time - 1., 360.*period_fac) //5 + 1.  
        data = data.groupby('xofyear').mean(('time'))
        u_mean, onset_pentad = scsm_onset(data.ucomp)
        onset.append(onset_pentad)
        data.close()
    return onset


load_onsets = np.load('/scratch/rg419/python_scripts/python_bin_updates/physics_updates/climatology/isca_era_onsets.npy')
onset_isca = np.asarray(load_onsets[0])
onset_era = np.asarray(load_onsets[1])

plt.hist(onset_isca)
plt.hist(onset_era, alpha=0.5)
plt.xlabel('Onset pentad')
plt.ylabel('Frequency')
plt.savefig(plot_dir+'isca_vs_era.pdf', format='pdf')
plt.close()
    

runs = ['control_qflux_0.500', 'control_qflux_0.750', 'control_qflux_1.250', 'control_qflux_1.500', 'control_qflux_1.750', 'control_qflux_2.000']
plt.hist(onset_isca)
for run in runs:
    print(run)
    onset = get_onsets(run)
    onset = [x for x in onset if x != None]
    np.save('onset_' + run, np.array(onset))
    if len(onset)!=0:
        plt.hist(onset, alpha=0.5)
        
plt.xlabel('Onset pentad')
plt.ylabel('Frequency')
plt.savefig(plot_dir+'isca_continents.pdf', format='pdf')
plt.close()



runs = ['rt_0.500', 'rt_0.750', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']

onset_1 = get_onsets('sn_1.000')
np.save('onset_sn_1.000', np.array(onset))
plt.hist(onset_1)
for run in runs:
    print(run)
    onset = get_onsets(run)
    onset = [x for x in onset if x != None]
    np.save('onset_' + run, np.array(onset))
    if len(onset)!=0:
        plt.hist(onset, alpha=0.5)

plt.xlabel('Onset pentad')
plt.ylabel('Frequency')
plt.savefig(plot_dir+'isca_rots.pdf', format='pdf')
plt.close()
    