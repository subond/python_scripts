"""Scripts for joining files etc."""

import os
import subprocess


#combine files with cdos
def join_files(var):
    
    nc_file_out = '/scratch/rg419/obs_and_reanalysis/era_' + var + '.nc'
    
    names = ['/scratch/rg419/obs_and_reanalysis/era_' + var + '_1979_1988.nc',
             '/scratch/rg419/obs_and_reanalysis/era_' + var + '_1989_1998.nc',
             '/scratch/rg419/obs_and_reanalysis/era_' + var + '_1999_2008.nc',
             '/scratch/rg419/obs_and_reanalysis/era_' + var + '_2009_2017.nc']
        
    nc_file_string = ' '.join(names)
        
    if not os.path.isfile(nc_file_out):
        subprocess.call('cdo -b F64 mergetime ' + nc_file_string + ' ' + nc_file_out, shell=True)
    

#Example
if __name__ == "__main__":
    join_files('w')

