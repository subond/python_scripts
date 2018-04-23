"""Scripts for joining files etc."""

import os
import subprocess
import sh


#combine files with cdos
def get_lev(var, lev):
    
    out_dir = '/scratch/rg419/obs_and_reanalysis/sep_levs_'+var+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(out_dir)
    
    nc_file_out = ['/scratch/rg419/obs_and_reanalysis/sep_levs_'+var+'/era_' + var + '_' + str(lev) + '_1979_1988.nc',
             '/scratch/rg419/obs_and_reanalysis/sep_levs_'+var+'/era_' + var + '_' + str(lev) + '_1989_1998.nc',
             '/scratch/rg419/obs_and_reanalysis/sep_levs_'+var+'/era_' + var + '_' + str(lev) + '_1999_2008.nc',
             '/scratch/rg419/obs_and_reanalysis/sep_levs_'+var+'/era_' + var + '_' + str(lev) + '_2009_2017.nc']
             
    names = ['/scratch/rg419/obs_and_reanalysis/original_era_data/era_' + var + '_1979_1988.nc',
             '/scratch/rg419/obs_and_reanalysis/original_era_data/era_' + var + '_1989_1998.nc',
             '/scratch/rg419/obs_and_reanalysis/original_era_data/era_' + var + '_1999_2008.nc',
             '/scratch/rg419/obs_and_reanalysis/original_era_data/era_' + var + '_2009_2017.nc']
    
    i=0
    for name in names:
        if not os.path.isfile(nc_file_out[i]):
            subprocess.call('cdo -sellevel,' + str(lev) + ' ' + name + ' ' + nc_file_out[i], shell=True)
        i=i+1


def join_files(var, lev):
    
    nc_file_out = '/scratch/rg419/obs_and_reanalysis/sep_levs_' + var + '/era_' + var + '_' + str(lev) + '.nc'
    
    names = ['/scratch/rg419/obs_and_reanalysis/sep_levs_'+var+'/era_' + var + '_' + str(lev) + '_1979_1988.nc',
             '/scratch/rg419/obs_and_reanalysis/sep_levs_'+var+'/era_' + var + '_' + str(lev) + '_1989_1998.nc',
             '/scratch/rg419/obs_and_reanalysis/sep_levs_'+var+'/era_' + var + '_' + str(lev) + '_1999_2008.nc',
             '/scratch/rg419/obs_and_reanalysis/sep_levs_'+var+'/era_' + var + '_' + str(lev) + '_2009_2017.nc']
        
    nc_file_string = ' '.join(names)
        
    if not os.path.isfile(nc_file_out):
        subprocess.call('cdo mergetime ' + nc_file_string + ' ' + nc_file_out, shell=True)
         

#Example
if __name__ == "__main__":
    for var in ['u']:
        for lev in range(350,1001,50):
            get_lev(var, lev)
            join_files(var, lev)
            

