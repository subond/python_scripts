from netCDF4 import Dataset
import numpy as np


def read_data( base_dir,exp_name,start_file,nfiles,avg_or_daily,thd_vars_read,twd_vars_read,ntime,ntime_tot,nlevs,nlats,nlons, model):
	"Function for reading in all data listed in thd_vars_read and twd_vars_read."
	nvars_thd=len(thd_vars_read)
	nvars_twd=len(twd_vars_read)
	
	#Create thd (3D) and twd (2D) data stores.
	thd_data=np.zeros((nvars_thd,ntime_tot,nlevs,nlats,nlons))
	twd_data=np.zeros((nvars_twd,ntime_tot,nlats,nlons))

	time_arr=np.zeros((ntime_tot))
	
	#For every file in the range, open it, read the 3D and 2D data, then close.
	for n in range(nfiles):
		nc_file_2d = base_dir+'/'+exp_name+'/run'+str(n+start_file)+'/atmos_'+avg_or_daily+'.nc'
		nc_file_3d = base_dir+'/'+exp_name+'/run'+str(n+start_file)+'/atmos_'+avg_or_daily+'_interp_model_lev.nc'
		print nc_file_2d
		fh_2d = Dataset(nc_file_2d, mode='r')
		fh_3d = Dataset(nc_file_3d, mode='r')
		for var_name in thd_vars_read:
			var_number=thd_vars_read.index(var_name)
			var_temp=fh_3d.variables[var_name][:]
# 			print thd_data.shape, var_temp.shape
# 			print thd_data[var_number,(n*ntime):((n+1)*ntime),:,:,:].shape
			thd_data[var_number,(n*ntime):((n+1)*ntime),:,:,:]=var_temp
			del var_temp
		for var_name in twd_vars_read:
			if var_name!='precip' or (model=='fms06' or model=='mima'):
				var_number=twd_vars_read.index(var_name)
				var_temp=fh_2d.variables[var_name][:]
				twd_data[var_number,(n*ntime):((n+1)*ntime),:,:]=var_temp
				del var_temp
			else:
				var_temp_1=fh_2d.variables['condensation_rain'][:]
				var_temp_2=fh_2d.variables['convection_rain'][:]
				var_temp=(var_temp_1+var_temp_2)*86400. #s 86400 is to change units of precip into mmday^-1. Because 1kgm^-2s^-1 = 86400 mm day^-1. 
				var_number=twd_vars_read.index(var_name)
				twd_data[var_number,(n*ntime):((n+1)*ntime),:,:]=var_temp				
				del var_temp,var_temp_1,var_temp_2

		time_arr[(n*ntime):((n+1)*ntime)] = fh_2d.variables['time'][:]
		fh_2d.close()
		fh_3d.close()

	idx_nan = (thd_data >= 1.0e+20)
	thd_data[idx_nan] = np.nan	
	return thd_data, twd_data, time_arr
	
	
def init( base_dir,exp_name,start_file,end_file,avg_or_daily,thd_vars_read,twd_vars_read ):
	"Uses the first nc file to read longitudes, lats etc."
	# File to initialise from
	nc_file_init = base_dir+'/'+exp_name+'/run'+str(start_file)+'/atmos_'+avg_or_daily+'.nc'
	nc_file_init_3d = base_dir+'/'+exp_name+'/run'+str(start_file)+'/atmos_'+avg_or_daily+'_interp_model_lev.nc'
# 	print nc_file_init, 'init directory'
	fh_init = Dataset(nc_file_init, mode='r')
	fh_init_3d = Dataset(nc_file_init_3d, mode='r')
	
	# Initialise variables that will be used throughout
	
	lons = fh_init.variables['lon'][:]
	lats = fh_init.variables['lat'][:]
	pfull = fh_init_3d.variables['pfull'][:]
#	time_bounds_init = fh_init.variables['time_bounds'][:] #This is specific and needs deleting once ntime has been found.
	time_init = fh_init.variables['time'][:] #This is specific and needs deleting once ntime has been found.

	
	#Initialise nlons, nlats, etc
	nlons=len(lons)
	nlats=len(lats)
	nlevs=len(pfull)
	ntime=len(time_init)
	nfiles=(end_file-start_file)+1
	ntime_tot=ntime*nfiles
	nvars_thd=len(thd_vars_read)
	nvars_twd=len(twd_vars_read)
	
# 	print time_init.shape
# 	print ntime
# 	print time_init

	
	
		
	# Delete any unwanted variables
	del time_init
	fh_init.close()
	fh_init_3d.close()
	

	
	return	lons,lats,pfull,nlons,nlats,nlevs,ntime,nfiles,ntime_tot,nvars_thd,nvars_twd

def read_land( base_dir,exp_name,nlats,nlons):
	"Function for reading in land mask."
	
	#Create thd (3D) and twd (2D) data stores.
	land_array=np.zeros((nlats,nlons))
	
	#For every file in the range, open it, read the 3D and 2D data, then close.
	nc_file = base_dir+'/exp/'+exp_name+'input/land_world_mountains.nc'
	print nc_file
	fh = Dataset(nc_file, mode='r')
	land_array=fh.variables['land_mask'][:]
	fh.close()

	return land_array

