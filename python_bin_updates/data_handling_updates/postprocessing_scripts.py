"""Scripts for joining files etc."""

import os
import subprocess

     
#combine files with cdos
def join_files(run, months, filename, filename_out = ' '):
    
    if filename_out == ' ':
        filename_out = filename + '_together'
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/' + filename + '.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    nc_file_string = ' '.join(names)
    
    nc_file_out = '/scratch/rg419/Data_moist/' + run + '/' + filename_out + '.nc'
    
    if not os.path.isfile(nc_file_out):
        subprocess.call('cdo mergetime ' + nc_file_string + ' ' + nc_file_out, shell=True)
    

#Example
if __name__ == "__main__":
    join_files('full_qflux', [121,481], 'tom_daily')

