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
    
    psi_mean = psi_mean.rolling(time=90).mean()
    
    sys.stderr.close()                # ordinary file object
    sys.stderr = temp
    
    return psi_mean

    
#sf_spinup('ss_eq_20', [[1,61]], filenames=['plev_monthly'])
#sf_spinup('ss_eq_20', [[1,61],[61,121]], filenames=['plev_monthly','plev_pentad'])
psi_2p5 = sf_spinup('ss_eq_2.5', [[61,121]], filenames=['plev_pentad'])
psi_5 = sf_spinup('ss_eq_5', [[61,121]], filenames=['plev_pentad'])
psi_10 = sf_spinup('ss_eq_10', [[61,121]], filenames=['plev_pentad'])
psi_15 = sf_spinup('ss_eq_15', [[61,121]], filenames=['plev_pentad'])
psi_20 = sf_spinup('ss_eq_20', [[61,121]], filenames=['plev_pentad'])

psi_2p5.plot.line()
psi_5.plot.line()
psi_10.plot.line()
psi_15.plot.line()
psi_20.plot.line()
plt.ylim([60,80])
plt.yticks(range(60,81,2))
plt.grid(True,linestyle=':')
plt.show()