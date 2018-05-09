"""
Check equilibration time of model using the vertically averaged overturning circulation strength. 2/05/2018
1) Check if model is spun up by looking at first 5 years of data
2) Check model response time to a sudden shift in insolation using second 5 years

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
from hadley_cell import mass_streamfunction
from data_handling_updates import cell_area#, rolling_mean
import sys
import statsmodels.api as sm


def sf_spinup(run, months_list, filenames=['plev_pentad']):
    
    # Function to open files for a specfied month range and filename.
    def open_files(run, months, filename):
        name_temp = '/scratch/rg419/Data_moist/' + run + '/run%04d/'+filename+'.nc'
        names = [name_temp % m for m in range( months[0], months[1])  ]
        data = xr.open_mfdataset( names, decode_times=False, chunks={'time': 30})
        # Reduce dataset so that only vcomp is retained
        data = xr.Dataset({'vcomp': data.vcomp}, coords=data.vcomp.coords) 
        #data.coords['month'] = data.time//30 + 1     
        #data = data.groupby('month').mean(('time'))
        return data
    
    arrays = []
    i=0
    for filename in filenames:
        data = open_files(run, months_list[i], filename)
        arrays.append(data)
        i=i+1
    
    data = xr.concat(arrays, dim='time')
        
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    
    area = cell_area(42, '/scratch/rg419/Isca/')
    area_xr = xr.DataArray(area, [('lat', data.lat ), ('lon', data.lon)])
    area_xr = area_xr.mean('lon')
    psi_mean = ((psi*area_xr).sum(('lat'))/area_xr.sum(('lat'))).mean('pfull')
    psi_mean = psi_mean*-1.
    
    # The rolling function will bring up alot of true_divide errors, but still works.
    # Can silence these by writing them to a file (although maybe this is slower)
    temp = sys.stderr                 # store original stdout object for later
    sys.stderr = open('log.txt', 'w') # redirect all prints to this log file
    
    #psi_mean = psi_mean.rolling(time=60).mean()
    #adj_time_70 = np.min(psi_mean.time.where(psi_mean >= 70.,drop=True).values)
    adj_time = np.min(psi_mean.time.where(psi_mean >= 74.,drop=True).values)
    #adj_time_75 = np.min(psi_mean.time.where(psi_mean >= 75.,drop=True).values)
    
    sys.stderr.close()                # ordinary file object
    sys.stderr = temp
    
    
    return psi_mean, adj_time

mlds = np.array([2.5,5.,10.,15.,20.])
adj_time_all = np.array([2292.5,2467.5,2642.5,2887.5,3192.5])-1800.

    
A = np.array([ mlds, np.ones(mlds.shape) ])

model = sm.OLS(adj_time_all, A.T)
result=model.fit()
consts = result.params
std_err = result.bse
        
print('=== Coeffs ===')
print(consts[0], consts[1])
print('=== Std Errs ===')
print(2*std_err[0], 2*std_err[1])


# Set plotting directory
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

plt.plot(mlds, adj_time_all, 'xk', mew=2, ms=10)
plt.plot(mlds, mlds * consts[0] + consts[1],'k')
plt.savefig(plot_dir + 'adj_time.pdf', format='pdf')
plt.close()

#sf_spinup('ss_eq_20', [[1,61]], filenames=['plev_monthly'])
psi_sst, adj_time = sf_spinup('ss_eq_sst', [[1,61],[61,121]], filenames=['plev_monthly','plev_pentad'])
#psi_sst, adj_time = sf_spinup('ss_eq_sst', [[61,121]], filenames=['plev_pentad'])
#psi_2p5, adj_time = sf_spinup('ss_eq_2.5', [[61,121]], filenames=['plev_pentad'])
#psi_5, adj_time = sf_spinup('ss_eq_5', [[61,121]], filenames=['plev_pentad'])
#psi_10, adj_time = sf_spinup('ss_eq_10', [[61,121]], filenames=['plev_pentad'])
#psi_15, adj_time = sf_spinup('ss_eq_15', [[61,121]], filenames=['plev_pentad'])
#psi_20, adj_time = sf_spinup('ss_eq_20', [[61,121]], filenames=['plev_pentad'])

psi_sst.plot.line()
#psi_2p5.plot.line()
#psi_5.plot.line()
#psi_10.plot.line()
#psi_15.plot.line()
#psi_20.plot.line()
plt.xlim([1700,1900])
#plt.ylim([60,80])
#plt.yticks(range(60,81,2))
plt.grid(True,linestyle=':')
plt.savefig(plot_dir + 'psi_sst_short.pdf', format='pdf')
plt.close()