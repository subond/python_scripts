import numpy as np
import os
from isca import IscaCodeBase, Experiment, GFDL_BASE
from isca.util import interpolate_output
import subprocess
from data_handling_updates import time_means
#from data_handling import join_files
 
#cb = IscaCodeBase.from_directory(GFDL_BASE)
    
#exp = Experiment('ss_91.000')
#for i in range(61, 248):
#     exp.runinterp(i,'atmos_pentad.nc','plev_pentad.nc', p_levs='EVEN')




for run in ["rt_0.500_5","rt_0.500_15", "rt_1.500_5","rt_1.500_15", "rt_1.750_5","rt_1.750_15","rt_2.000_5",]: 
    print(run)
    for i in range(121,481):
        #print(i)
        try:
            infile = '/scratch/rg419/Data_moist/' + run + '/run%04d/atmos_pentad.nc' % i
            outfile = '/scratch/rg419/Data_moist/' + run + '/run%04d/plev_pentad.nc' % i
            interpolate_output(infile, outfile, p_levs='EVEN', var_names=['slp', 'height'])
        
            infile = '/scratch/rg419/Data_moist/' + run + '/run%04d/atmos_daily.nc' % i
            outfile = '/scratch/rg419/Data_moist/' + run + '/run%04d/plev_daily.nc' % i
            interpolate_output(infile, outfile, p_levs='EVEN')
        
        except:
            print(i)
    test=time_means(run, [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)

#test=time_means('rt_0.750_5', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)
#test=time_means('rt_1.250_5', [121,481], filename='plev_pentad', timeav='pentad', write_netcdf=True)


        #infile = '/scratch/rg419/Data_moist/' + run + '/run%04d/atmos_daily_mean.nc' % i
        #outfile = '/scratch/rg419/Data_moist/' + run + '/run%04d/plev_daily_mean.nc' % i
        #interpolate_output(infile, outfile, p_levs='EVEN')
        
        

#for run in ["sine_sst_10m_ss_100", "sine_sst_10m_ss_110", "sine_sst_10m_ss_120", "sine_sst_10m_ss_130",
#            "sine_sst_10m_ss_140", "sine_sst_10m_ss_150", "sine_sst_10m_ss_160", "sine_sst_10m_ss_170", "sine_sst_10m_ss_180"]: 
#    exp = Experiment(run, codebase=cb)
#    for i in range(61, 181):
#        runinterp(exp, i,'atmos_pentad.nc','plev_pentad.nc', p_levs='EVEN', var_names='-a slp height')
#        runinterp(exp, i,'atmos_daily.nc','plev_daily.nc', p_levs='EVEN')
