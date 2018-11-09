'''01/10/2018 Load up onset dates for ERA-Interim, JRA-55, and NCEP/NCAR and compare between datasets, and compare with dates from published studies'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm
from data_handling_updates import gradients as gr

plot_dir = '/scratch/rg419/plots/onset_variability/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

# Load in onset dates
onsets_scsm_ncep = np.load('/scratch/rg419/python_scripts/onset_variability/ncep_onsets_scsm.npy')
onsets_scsm_ncep = xr.DataArray(onsets_scsm_ncep, coords={'year': ('year', range(1948,2017))}, dims=['year'])
onsets_bob_ncep = np.load('/scratch/rg419/python_scripts/onset_variability/ncep_onsets_bob_days.npy')
onsets_bob_ncep = xr.DataArray(onsets_bob_ncep, coords={'year': ('year', range(1948,2017))}, dims=['year'])
onsets_ism_ncep = np.load('/scratch/rg419/python_scripts/onset_variability/ncep_onsets_ism_days.npy')
onsets_ism_ncep = xr.DataArray(onsets_ism_ncep, coords={'year': ('year', range(1948,2017))}, dims=['year'])

onsets_scsm_jra = np.load('/scratch/rg419/python_scripts/onset_variability/jra_onsets_scsm.npy')
onsets_scsm_jra = xr.DataArray(onsets_scsm_jra, coords={'year': ('year', range(1958,2017))}, dims=['year'])
onsets_bob_jra = np.load('/scratch/rg419/python_scripts/onset_variability/jra_onsets_bob_days.npy')
onsets_bob_jra = xr.DataArray(onsets_bob_jra, coords={'year': ('year', range(1958,2017))}, dims=['year'])
onsets_ism_jra = np.load('/scratch/rg419/python_scripts/onset_variability/jra_onsets_ism_days.npy')
onsets_ism_jra = xr.DataArray(onsets_ism_jra, coords={'year': ('year', range(1958,2017))}, dims=['year'])

onsets_scsm_era = np.load('/scratch/rg419/python_scripts/onset_variability/era_onsets_scsm.npy')
onsets_scsm_era = xr.DataArray(onsets_scsm_era, coords={'year': ('year', range(1979,2017))}, dims=['year'])
onsets_bob_era = np.load('/scratch/rg419/python_scripts/onset_variability/era_onsets_bob_days.npy')
onsets_bob_era = xr.DataArray(onsets_bob_era, coords={'year': ('year', range(1979,2017))}, dims=['year'])
onsets_ism_era = np.load('/scratch/rg419/python_scripts/onset_variability/era_onsets_ism_days.npy')
onsets_ism_era = xr.DataArray(onsets_ism_era, coords={'year': ('year', range(1979,2017))}, dims=['year'])

#SCS onset Wang et al 2004
wang_onset_scsm = np.array([26,30,26,25,27,26,31,29,31,32,29,30,30,27,28,30,28,29,25,29,34,29,32,25,26,33,29,
                            31,26,28,29,27,27,31,31,31,29,30,27,32,29,28,28,32,28,32,25,27,26,28,29,30,26,26])

wang_onset_scsm = xr.DataArray(wang_onset_scsm, coords={'year': ('year', range(1948,2002))}, dims=['year'])

#ISM onset Wang et al 2009
wang_onset_ism = np.array([36,25,37,32,23,38,31,31,19,31,39,27,15,19,18,34,37,36,35,40,39,24,27,25,46,36,27,
                           31,29,38,29,44,33,33,32,45,29,23,37,33,32,21,17,33,39,29,30,39,33,46,37,19,15,24,
                           35,38,19,37,23,42]) + 120

wang_onset_ism = xr.DataArray(wang_onset_ism, coords={'year': ('year', range(1948,2008))}, dims=['year'])

# BOB onset Mao and Wu 2007
maowu_onset_bob = np.array([123,108,132,122,120,145,123,125,118,115,118,130,120,115,120,122,106,114,117,131,134,129,
                            128,129,112,139,106,108,128,122,121,111,131,115,135,136,121,129,114,136,135,101,104,117])

maowu_onset_bob = xr.DataArray(maowu_onset_bob, coords={'year': ('year', range(1958,2002))}, dims=['year'])




# Plot onsets 
plt.plot(onsets_scsm_ncep.year, onsets_scsm_ncep,'o-', color='C0')
plt.plot(onsets_scsm_jra.year, onsets_scsm_jra,'x-.', color='C1')
plt.plot(onsets_scsm_era.year, onsets_scsm_era,'+--', color='C2')
plt.plot(wang_onset_scsm.year, wang_onset_scsm,'o:', markerfacecolor='none', color='k')
plt.xlabel('Year')
plt.ylabel('Onset pentad')
plt.yticks(range(22,40,2))
plt.grid(True,linestyle=':')
plt.legend(['NCEP/NCAR','JRA-55','ERA-Interim','Wang et al. 2004'])
plt.title('South China Sea Monsoon Onset')
plt.savefig(plot_dir + 'scs_onsets_reanalysis.pdf', format='pdf')
plt.close()



plt.plot(onsets_bob_ncep.year, onsets_bob_ncep,'o-', color='C0')
plt.plot(onsets_bob_jra.year, onsets_bob_jra,'x-.', color='C1')
plt.plot(onsets_bob_era.year, onsets_bob_era,'+--', color='C2')
plt.plot(maowu_onset_bob.year, maowu_onset_bob,'o:', markerfacecolor='none', color='k')
plt.xlabel('Year')
plt.ylabel('Onset day')
plt.yticks(range(90,156,5))
plt.grid(True,linestyle=':')
plt.legend(['NCEP/NCAR','JRA-55','ERA-Interim','Wang et al. 2004'])
plt.title('Bay of Bengal Monsoon Onset')
plt.savefig(plot_dir + 'bob_onsets_reanalysis.pdf', format='pdf')
plt.close()



plt.plot(onsets_ism_ncep.year, onsets_ism_ncep,'o-', color='C0')
plt.plot(onsets_ism_jra.year, onsets_ism_jra,'x-.', color='C1')
plt.plot(onsets_ism_era.year, onsets_ism_era,'+--', color='C2')
plt.plot(wang_onset_ism.year, wang_onset_ism,'o:', markerfacecolor='none', color='k')
plt.xlabel('Year')
plt.ylabel('Onset day')
plt.yticks(range(130,175,5))
plt.grid(True,linestyle=':')
plt.legend(['NCEP/NCAR','JRA-55','ERA-Interim','Wang et al. 2004'])
plt.title('Indian Summer Monsoon Onset')
plt.savefig(plot_dir + 'ism_onsets_reanalysis.pdf', format='pdf')
plt.close()


