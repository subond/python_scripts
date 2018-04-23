"""Script to take model output and produce a climatology. 
   23/11/2017 modified to allow files in hours to be converted to days using day_fac
   NB It's good to check that rundata.xofyear contains the values you'd expect"""

import numpy as np
import xarray as xr
import scipy.interpolate as spint

def interp_datafiles(run, months, name_in='daily', name_out='plev_daily'):
    
    offline_interp = ['div', 'height', 'omega', 'temp', 'teq', 'ucomp', 'vcomp', 'vor']
    
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+name_in+'.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    outfile_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+name_out+'.nc'
    outfiles = [outfile_temp % m for m in range( months[0], months[1])  ]
    
    for i in range(len(names)):
        data = xr.open_dataset( names[i], decode_times=False)
        pfull_new = np.arange(1000., 0., -50.)        
        data_interp = xr.Dataset({'ps': (['time', 'lat', 'lon'], data.ps),
                     'h_trop':  (['time', 'lat', 'lon'], data.h_trop)},
                     coords={'time': ('time', data.time),
                             'pfull': ('pfull', pfull_new),
                               'lat': ('lat', data.lat),
                               'lon': ('lon', data.lon),
                               'latb': ('latb', data.latb),
                               'lonb': ('lonb', data.lonb)})
        for data_var in offline_interp:
            print data_var
            f = spint.interp1d(data.pfull, data[data_var], axis=1, fill_value='extrapolate', kind='quadratic')
            data_interp[data_var] = xr.DataArray(f(pfull_new), coords=[data.time, pfull_new, data.lat, data.lon], dims=['time', 'pfull', 'lat', 'lon'])
        
        data = data_interp
        
        filename = outfiles[i]
        data.to_netcdf(filename)
    

    

if __name__ == "__main__":
    
    test=interp_datafiles('dry_zs', [9,21])