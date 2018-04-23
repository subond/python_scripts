import numpy as np
import os
from gfdl import Experiment
import subprocess
from data_handling import join_files
 
     
#exp = Experiment('ss_91.000')
#for i in range(61, 248):
#     exp.runinterp(i,'atmos_pentad.nc','plev_pentad.nc', p_levs='EVEN')

for run in ["zs_sst", "ap10_co2", "ap10_qflux"]: 
    exp = Experiment(run)
    for i in range(121, 481):
        exp.runinterp(i,'atmos_pentad.nc','plev_pentad.nc', p_levs='EVEN', var_names='-a slp height')
        exp.runinterp(i,'atmos_daily.nc','plev_daily.nc', p_levs='EVEN')
        
#for run in ["qflux_0_200_small"]: 
#    exp = Experiment(run)
#    for i in range(121, 241):
#        exp.runinterp(i,'atmos_pentad.nc','plev_pentad.nc', p_levs='EVEN')
#        exp.runinterp(i,'atmos_daily.nc','plev_daily.nc', p_levs='EVEN')

#p_levs=' -p "10 30 100 300 500 700 1000 3000 5000 7000 10000 15000 20000 25000 30000 40000 50000 60000 70000 75000 80000 85000 90000 95000 100000" '

#exp = Experiment('full_qflux')
#for i in range(121, 481):
#     exp.runinterp(i,'atmos_pentad.nc','tom_ptd.nc', p_levs=p_levs, var_names='height hght', extra_flags='-x')
     
#join_files('full_qflux', [121,481], 'atmos_daily')