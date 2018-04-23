# Load up experiment files and plot spinup
# Inputs:
# run: name of run
# field: specify field to plot
# months_list: list of 2 element lists of start and end months
# filenames: list of names of individual files for a given month range
# plevs: Used for 3D fields. List of: top pressure level, bottom pressure level, name to refer to this range as. Default is to do all levels

# NB Doesn't work for 6 hourly data yet

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import cell_area

def spinup_fn(run, field, months_list, filenames=['atmos_pentad'], plevs=[0.,2000.,'all']):
    
    # Function to open files for a specfied month range and filename.
    # Takes annual means
    def open_files(run, months, filename):
        name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
        names = [name_temp % m for m in range( months[0], months[1])  ]
        #read data into xarray 
        data = xr.open_mfdataset( names, decode_times=False, chunks={'time': 30})
        data.coords['year'] = data.time//360 + 1        
        field_yr = data[field].groupby('year').mean(('time'))
        
        return field_yr, data
    
    # Combine data from files with different names (eg. atmos_monthly and atmos_pentad) into one time series
    arrays = []
    i=0
    for filename in filenames:
        field_yr,data = open_files(run, months_list[i], filename)
        arrays.append(field_yr)
        i=i+1
    
    field_yr = xr.concat(arrays, dim='year')
    
    # Check if data is 3D and if so integrate over specfied levels
    try:
        p_levs = data.pfull[ (data.pfull >= plevs[0]) & (data.pfull <= plevs[1])]
        dp = xr.DataArray( np.diff(data.phalf), [('pfull', field_yr.pfull) ]) * 100.
        field_yr = (field_yr*dp).sel(pfull=p_levs).sum('pfull')/9.8
        print '3D field, vertical integral taken'
        three_d = True
    
    except:
        print '2D field'
        three_d = False
    
    # Calculate cell areas and take area mean
    area = cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')
    area_xr = xr.DataArray(area, [('lat', data.lat ), ('lon', data.lon)])
    field_av = (field_yr*area_xr).sum(('lat','lon'))/area_xr.sum(('lat','lon'))
    
    # Plot up result and save
    field_av.plot()
    plt.xlabel('Year')
    plt.ylabel(field)
    if three_d:
        plotname = '/scratch/rg419/plots/spinup/'+ field + '_' + str(plevs[2]) + '_spinup_'+ run +'.png'
    else:
        plotname = '/scratch/rg419/plots/spinup/'+ field + '_spinup_'+ run +'.png'
    plt.savefig(plotname)
    plt.close()
    

    return field_av



if __name__ == "__main__":

	#t_av = spinup_fn('ss_90.000', 't_surf', months_list= [[1,61],[61,181]], filenames=['atmos_monthly', 'atmos_pentad'])
	#t_av = spinup_fn('ss_95.000', 't_surf', months_list= [[1,61],[61,181]], filenames=['atmos_monthly', 'atmos_pentad'])
	#t_av = spinup_fn('ss_100.000', 't_surf', months_list= [[1,61],[61,181]], filenames=['atmos_monthly', 'atmos_pentad'])
	#t_av = spinup_fn('ss_100.000', 'sphum', months_list= [[1,61],[61,181]], filenames=['atmos_monthly', 'atmos_pentad'])
	#spinup_fn('qflux_25_300', 't_surf', months_list= [[1,121],[121,218]], filenames=['atmos_monthly', 'atmos_pentad'])
	spinup_fn('qflux_25_300', 'sphum', months_list= [[1,121],[121,218]], filenames=['atmos_monthly', 'atmos_pentad'])
	spinup_fn('qflux_25_300', 'sphum', months_list= [[1,121],[121,218]], filenames=['atmos_monthly', 'atmos_pentad'], plevs=[0.,100.,'strat'])

