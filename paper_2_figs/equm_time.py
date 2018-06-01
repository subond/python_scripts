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
    adj_time = np.min(data.p_cent.time.where(data.p_cent >= 15.,drop=True).values)    
    
    return data.p_cent, adj_time


# Set plotting directory
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


#sf_spinup('ss_eq_20', [[1,61]], filenames=['plev_monthly'])
#psi_sst, adj_time = sf_spinup('ss_eq_sst', [[1,61],[61,121]], filenames=['plev_monthly','plev_pentad'])
pcent_sst, adj_time_sst = sf_spinup('ss_eq_sst', [[61,121]], filenames=['plev_pentad'])
print('sst')
pcent_sst_zs, adj_time_sst_zs = sf_spinup('ss_eq_sst_zs', [[61,121]], filenames=['plev_pentad'])
print('sst_zs')
pcent_2p5, adj_time_2p5 = sf_spinup('ss_eq_2.5', [[61,121]], filenames=['plev_pentad'])
print('2.5')
pcent_5, adj_time_5 = sf_spinup('ss_eq_5', [[61,121]], filenames=['plev_pentad'])
print('5')
pcent_10, adj_time_10 = sf_spinup('ss_eq_10', [[61,121]], filenames=['plev_pentad'])
print('10')
pcent_15, adj_time_15 = sf_spinup('ss_eq_15', [[61,121]], filenames=['plev_pentad'])
print('15')
pcent_20, adj_time_20 = sf_spinup('ss_eq_20', [[61,121]], filenames=['plev_pentad'])
print('20')


# Set figure parameters
rcParams['figure.figsize'] = 5, 6
rcParams['font.size'] = 16

fig, (ax1, ax2) = plt.subplots(2)

sst_zs, = ax1.plot(pcent_sst_zs.time-1800., pcent_sst_zs,'y')
sst, = ax1.plot(pcent_sst.time-1800., pcent_sst,'m')
m25, = ax1.plot(pcent_2p5.time-1800., pcent_2p5, 'b', label='2.5')
m5, = ax1.plot(pcent_5.time-1800., pcent_5, 'g', label='5.')
m10, = ax1.plot(pcent_10.time-1800., pcent_10, 'k', label='10.')
m15, = ax1.plot(pcent_15.time-1800., pcent_15, 'r', label='15.')
m20, = ax1.plot(pcent_20.time-1800., pcent_20, 'c', label='20.')
ax1.set_ylabel('Precipitation centroid')
ax1.set_xlabel('Time, days')
ax1.set_xlim([0,1000])
ax1.set_ylim([-2.5,25])
ax1.set_yticks(range(0,26,5))
ax1.grid(True,linestyle=':')
legend = ax1.legend([m25,m5,m10,m15,m20,sst,sst_zs], ['2.5','5.','10.','15.','20.','SST','SST (ZS)'], loc='lower right', fontsize=8, title='MLD, m', ncol=2) #, ,#bbox_to_anchor=(1.05, 1),
#legend = ax1.legend([sst,sst_zs], ['SST','SST (ZS)'], loc='lower right', fontsize=8) #, ,#bbox_to_anchor=(1.05, 1),
legend.get_title().set_fontsize(8)

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

plt.savefig(plot_dir + 'eq_time.pdf', format='pdf')
plt.close()