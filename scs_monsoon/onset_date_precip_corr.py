""" 14/06/2018 Look at the correlation between the SCS onset date and the precipitation amount, plot as a map
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams

rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 14

plot_dir = '/scratch/rg419/plots/scs_monsoon/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


load_onsets = np.load('/scratch/rg419/python_scripts/python_bin_updates/physics_updates/climatology/isca_era_onsets.npy')
onset_era = np.asarray(load_onsets[1])

def get_mjjas_precip(years):
    precip_mjjas = []
    for year in years:
        name_temp = '/scratch/rg419/obs_and_reanalysis/GPCP_monthly/gpcp_cdr_v23rB1_y' + str(year) + '_m%02d.nc'
        names = [name_temp % m for m in range(5,10)]
        #read data into xarray 
        data = xr.open_mfdataset(names)
        precip_mjjas.append(data.precip.mean('time').values)
        lats = data.precip.latitude
        lons = data.precip.longitude
    return np.asarray(precip_mjjas), lats, lons

precip_mjjas, lats, lons = get_mjjas_precip(range(1979,2017))

precip_mjjas = xr.DataArray(precip_mjjas, coords=[range(1979,2017), lats.values, lons.values], dims=['year', 'lat', 'lon'])
onset_era = xr.DataArray(onset_era, coords=[range(1979,2017)], dims=['year'])

onset_mean = onset_era.mean('year')
onset_diffs = onset_era - onset_mean

precip_mean = precip_mjjas.mean('year')
precip_diffs = precip_mjjas - precip_mean

cov_precip_onset = (onset_diffs * precip_diffs).mean('year')
std_precip = np.sqrt((precip_diffs**2.).mean('year'))
std_onset = np.sqrt((onset_diffs**2.).mean('year'))

corr = cov_precip_onset/std_precip/std_onset

land_data = xr.open_dataset('/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc')

corr.plot(levels=[-0.4,-0.3,-0.2,0.2,0.3,0.4], add_labels=False)
land_data.lsm[0,:,:].plot.contour(levels=[0.,100.], colors='k', add_labels=False)
#plt.xlim([70,155])
#plt.ylim([-15,60])
plt.title('Correlation of precip and SCS onset pentad')
plt.xlabel('Latitude')
plt.ylabel('Longitude')
plt.savefig(plot_dir+'onset_date_precip_corr_fullmap.pdf', format='pdf')
