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
from data_handling_updates import cell_area, gradients as gr
import sys
import statsmodels.api as sm
from climatology import precip_centroid

def sf_spinup(run, months_list, filenames=['plev_pentad']):
    
    # Function to open files for a specfied month range and filename.
    def open_files(run, months, filename):
        name_temp = '/disca/share/rg419/Data_moist/' + run + '/run%04d/'+filename+'.nc'
        names = [name_temp % m for m in range( months[0], months[1])  ]
        data = xr.open_mfdataset( names, decode_times=False, chunks={'time': 30})
        # Reduce dataset so that only vcomp is retained
        data = xr.Dataset({'vcomp': data.vcomp, 'precipitation':data.precipitation}, coords=data.vcomp.coords) 
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
    
    precip_centroid(data)
        
    #psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    #psi /= 1.e9
    
    #area = cell_area(42, '/scratch/rg419/Isca/')
    #area_xr = xr.DataArray(area, [('lat', data.lat ), ('lon', data.lon)])
    #area_xr = area_xr.mean('lon')
    #psi_mean = ((psi*area_xr).sum(('lat'))/area_xr.sum(('lat'))).mean('pfull')
    #psi_mean = psi_mean*-1.
    
    # The rolling function will bring up alot of true_divide errors, but still works.
    # Can silence these by writing them to a file (although maybe this is slower)
    temp = sys.stderr                 # store original stdout object for later
    sys.stderr = open('log.txt', 'w') # redirect all prints to this log file
    
    #psi_mean = gr.ddt(psi_mean, timedir='time', secperunit=86400., cyclic=False)
    
    #psi_mean = psi_mean.rolling(time=6).mean()
    psi_mean = data.p_cent
    #adj_time_70 = np.min(psi_mean.time.where(psi_mean >= 70.,drop=True).values)
    adj_time = np.min(psi_mean.time.where(psi_mean >= 15.,drop=True).values)
    #adj_time_75 = np.min(psi_mean.time.where(psi_mean >= 75.,drop=True).values)
    #adj_time=1.
    sys.stderr.close()                # ordinary file object
    sys.stderr = temp
    
    
    return psi_mean, adj_time


# Set plotting directory
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


#sf_spinup('ss_eq_20', [[1,61]], filenames=['plev_monthly'])
#psi_sst, adj_time = sf_spinup('ss_eq_sst', [[1,61],[61,121]], filenames=['plev_monthly','plev_pentad'])
psi_sst, adj_time_sst = sf_spinup('ss_eq_sst', [[61,121]], filenames=['plev_pentad'])
print('sst')
psi_sst_zs, adj_time_sst_zs = sf_spinup('ss_eq_sst_zs', [[61,121]], filenames=['plev_pentad'])
print('sst_zs')
psi_2p5, adj_time_2p5 = sf_spinup('ss_eq_2.5', [[61,121]], filenames=['plev_pentad'])
print('2.5')
psi_5, adj_time_5 = sf_spinup('ss_eq_5', [[61,121]], filenames=['plev_pentad'])
print('5')
psi_10, adj_time_10 = sf_spinup('ss_eq_10', [[61,121]], filenames=['plev_pentad'])
print('10')
psi_15, adj_time_15 = sf_spinup('ss_eq_15', [[61,121]], filenames=['plev_pentad'])
print('15')
psi_20, adj_time_20 = sf_spinup('ss_eq_20', [[61,121]], filenames=['plev_pentad'])
print('20')


# Set figure parameters
rcParams['figure.figsize'] = 5, 6
rcParams['font.size'] = 16

fig, (ax1, ax2) = plt.subplots(2)

#ax1.plot(psi_sst.time-1800., psi_sst)
#ax1.plot(psi_sst_zs.time-1800., psi_sst_zs)
ax1.plot(psi_2p5.time-1800., psi_2p5)
ax1.plot(psi_5.time-1800., psi_5)
ax1.plot(psi_10.time-1800., psi_10)
ax1.plot(psi_15.time-1800., psi_15)
ax1.plot(psi_20.time-1800., psi_20)
ax1.set_ylabel('Precipitation centroid')
ax1.set_xlabel('Time, days')
ax1.set_xlim([0,1000])
ax1.set_ylim([-2.5,40])
ax1.set_yticks(range(0,21,5))
ax1.grid(True,linestyle=':')

mlds = np.array([2.5,5.,10.,15.,20.])
adj_time_all = np.array([adj_time_2p5, adj_time_5, adj_time_10, adj_time_15, adj_time_20]) - 1800.

A = np.array([ mlds, np.ones(mlds.shape) ])

model = sm.OLS(adj_time_all, A.T)
result=model.fit()
consts = result.params
std_err = result.bse
        
print('=== Coeffs ===')
print(consts[0], consts[1])
print('=== Std Errs ===')
print(2*std_err[0], 2*std_err[1])

plt.plot(mlds, adj_time_all, 'xk', mew=2, ms=10)
plt.plot(mlds, mlds * consts[0] + consts[1],'k')
ax2.set_ylabel('Response time, days')
ax2.set_xlabel('Mixed layer depth, m')
ax2.grid(True,linestyle=':')
ax2.set_yticks(range(200,801,200))

plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1, hspace=0.3)

plt.savefig(plot_dir + 'adj_time_30win.pdf', format='pdf')
plt.close()