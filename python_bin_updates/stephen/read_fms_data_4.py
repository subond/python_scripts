import sys
from netCDF4 import Dataset
import numpy as np
import nc_file_io_xarray as io
import plotting_data as splot
import analyse as ana
import thetaint_version4_2015_earth as spv
import matplotlib.pyplot as plt
import xarray as xarray
# model options are fms06, mima or fms13.
model='fms13'
land_present = True
topo_present = True

just_plot_and_read_pickle_data = False


if model=='fms13':	
	input_dir='/scratch/sit204/FMS2013/GFDLmoistModel/'
	base_dir='/scratch/sit204/Data_2013/'
	land_file='input/land_world_mountains.nc'
	base_exp_name='warmpool_exp_mk_1/'
	exp_name='warmpool_exp_mk1_1/'
#	exp_name_2='warmpool_cool_1/'
else:
	print 'Invalid model type.'
	sys.exit()

avg_or_daily='6hourly'
ptemp_calc=False
brunt_vas_calc=False
do_pv_calc=False
do_cfl_calc=True

start_file=13
end_file=116

brunt_vas_level_1 = 850. #Levels to measure brunt vasala between (in hPa)
brunt_vas_level_2 = 925. #Levels to measure brunt vasala between (in hPa)

data_file_name = 'python_data_store'

#**** end option setting ********


planet_params = ana.planet_params_set()
model_params = ana.model_params_set(input_dir)
outfile = base_dir + exp_name + data_file_name #Set the location of the output data file that will be saved / read from
nfiles=(end_file-start_file)+1

analysis_list={}
analysis_list['ptemp']=ptemp_calc
analysis_list['brunt_vas']=brunt_vas_calc
analysis_list['pv']=do_pv_calc
analysis_list['cfl']=do_cfl_calc

#***********Read data*************************
dataset, time_arr, size_list = io.read_data( base_dir,exp_name,start_file,end_file,avg_or_daily,topo_present)

land_array = io.read_land(input_dir,base_exp_name,land_present,size_list,land_file)
dataset['land'] = (('lat','lon'),land_array)

#s If there is a comparison to be made, try reading from second experiment name.

try:
	exp_name_2
except NameError:
	print 'running analysis on one experiment only'
else:
	dataset_2, time_arr2, size_list = io.read_data( base_dir,exp_name_2,start_file,end_file,avg_or_daily,topo_present)
	dataset_2['land'] = (('lat','lon'),land_array)

	if((time_arr2 == time_arr).to_index().all()):
		dataset_diff=dataset-dataset_2
		dataset_diff['land'] =(('lat','lon'),land_array)
		dataset_diff.attrs['exp_name']='diff_'+exp_name+'_minus_'+exp_name_2
#************End Read data********************


#***********Analyse data*************************
dataset = ana.run_analysis(analysis_list, dataset, model_params)
#***********End analyse data*********************

#***********Plotting data*************************
plot_list={}

#plot_list[len(plot_list.keys())]={'var_name':'ucomp','plot_type':'latp', 'av_type':('lon','time'), 'group_type':'years', 'y_axis_type':'logp'}
plot_list[len(plot_list.keys())]={'var_name':'cfl','plot_type':'latp', 'av_type':('lon','time'), 'group_type':'years', 'y_axis_type':'logp'}

plot_list[len(plot_list.keys())]={'var_name':'cfl','plot_type':'latlon','av_type':('time'), 'group_type':'months', 'plev':0.01}

#plot_list[len(plot_list.keys())]={'var_name':'ucomp','plot_type':'latlon','av_type':('time'), 'group_type':'months', 'plev':850.}


#plot_list[len(plot_list.keys())]={'var_name':'t_surf','plot_type':'latlon','av_type':('time'), 'group_type':'seasons'}

#plot_list[len(plot_list.keys())]={'var_name':'t_surf','plot_type':'zmlat','av_type':('lon','time'), 'group_type':'all_time'}

splot.run_plotting(dataset, plot_list)
#splot.run_plotting(thd_data_2,twd_data_2, plot_list)
#splot.run_plotting(thd_data_diff,twd_data_diff, plot_list)

#***********End plotting data*************************








