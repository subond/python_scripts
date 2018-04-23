import xarray as xarray
from netCDF4 import Dataset
import numpy as np
import calendar_calc as cal
import os
import sys
import pdb
import create_timeseries as cts

def read_data( base_dir,exp_name,start_file,end_file,avg_or_daily,topo_present):

	files_temp=[base_dir+'/'+exp_name+'/run%d' % m for m in range(start_file, end_file+1)]
	if(topo_present):
		extra='_interp_new_model_lev.nc'
	else:
		extra='.nc'

	thd_string = '/atmos_'+avg_or_daily+extra

	thd_files = [s + thd_string for s in files_temp]

	thd_files_exist=[os.path.isfile(s) for s in thd_files]

	print thd_files[0]

	if not(all(thd_files_exist)):
		print 'WARNING missing files', 	[thd_files[elem] for elem in range(len(thd_files_exist)) if not thd_files_exist[elem]]
		print 'EXITING BECAUSE OF MISSING FILES'
		sys.exit(0)

	size_list = init(thd_files[0])

	da_3d = xarray.open_mfdataset(	thd_files,
			decode_times=False,  # no calendar so tell netcdf lib
			# choose how data will be broken down into manageable
			# chunks.
#                        concat_dim='time',
			chunks={'time': size_list['ntime'],
					'lon': size_list['nlons']//4,
					'lat': size_list['nlats']//2})

	time_arr = da_3d.time
	date_arr = cal.day_number_to_date(time_arr)

	da_3d.coords['dayofyear_ax'] = (('dayofyear_ax'),np.unique(np.mod(np.floor(time_arr),360)))
	da_3d.coords['months_ax'] = (('months_ax'),np.unique(date_arr.month))
	da_3d.coords['seasons_ax'] = (('seasons_ax'),np.arange(4))
#	da_3d.coords['years_ax'] = (('years_ax'),date_arr.year)

	da_3d.coords['dayofyear'] = (('time'),np.mod(np.floor(time_arr),360))
	da_3d.coords['months'] = (('time'),date_arr.month)
	da_3d.coords['years'] = (('time'),date_arr.year)
	da_3d.coords['seasons'] = (('time'),np.mod((da_3d.time +30.) // 90,4))
	da_3d.coords['offset_seasons'] = (('time'),np.mod((da_3d.time -30.) // 90,4))

	da_3d.coords['all_time'] = (('time'),time_arr/time_arr)

	da_3d.coords['seq_months'] = (('time'),time_arr//30 + 1)
	da_3d.coords['seq_seasons'] = (('time'),(da_3d.time +30.) // 90)

	da_3d.coords['seq_seasons_ax'] = (('seq_seasons_ax'),np.mod(np.min(da_3d.seq_seasons.values),4)+np.arange(len(np.unique(da_3d.seq_seasons.values))))

	da_3d.attrs['exp_name']=exp_name
	try:
		da_3d['precipitation']
	except KeyError:
		try:
		   print 'aggregating rain'
		   da_3d['convection_rain']
		   da_3d['condensation_rain']
		except KeyError:
		   print 'no precip output present'
		else:		
                   da_3d['precip']=(('time','lat','lon'),da_3d['convection_rain']+da_3d['condensation_rain'])
		   print 'done aggregating rain'




	thd_data = da_3d

	return thd_data, time_arr, size_list

def init( nc_file_init):
	"Uses the first nc file to read longitudes, lats etc."
	fh_init = Dataset(nc_file_init, mode='r')
	
	# Initialise variables that will be used throughout
	
	lons = fh_init.variables['lon'][:]
	lats = fh_init.variables['lat'][:]
	pfull = fh_init.variables['pfull'][:]
	time_init = fh_init.variables['time'][:] #This is specific and needs deleting once ntime has been found.

	
	#Initialise nlons, nlats, etc
	nlons=len(lons)
	nlats=len(lats)
	nlevs=len(pfull)
	ntime=len(time_init)
	
	# Delete any unwanted variables
	del time_init
	fh_init.close()
	
	size_list={}
	size_list['nlons']=nlons
	size_list['nlats']=nlats
	size_list['nlevs']=nlevs
	size_list['ntime']=ntime
	
	return	size_list

def read_land( base_dir,exp_name,land_present,topo_present,size_list,land_file='input/land.nc'):
	"Function for reading in land mask."
	
	#Create thd (3D) and twd (2D) data stores.
	land_array=np.zeros((size_list['nlats'],size_list['nlons']))
	topo_array=np.zeros((size_list['nlats'],size_list['nlons']))
	
	if(land_present or topo_present):
		nc_file = base_dir+'/exp/'+exp_name+land_file
		print nc_file
		fh = Dataset(nc_file, mode='r')
		if land_present:
			land_array=fh.variables['land_mask'][:]
		if topo_present:
			topo_array=fh.variables['zsurf'][:]
		fh.close()

	return land_array, topo_array

def diff_time_adjust(dataset,dataset_2):

	ntime_1=len(dataset.time.to_index())
	ntime_2=len(dataset_2.time.to_index())

	if ntime_1!=ntime_2:
		print 'FATAL - total number of days in dataset_1 and dataset_2 are different. In dataset 1 and 2 there are ',ntime_1,ntime_2,' times, respectively.'
		sys.exit()

	time_arr_1=dataset.time.to_index()
	time_arr_2=dataset_2.time.to_index()

	doy_arr_1=dataset.dayofyear.to_index()
	doy_arr_2=dataset_2.dayofyear.to_index()

	if(all(time_arr_1!=time_arr_2)):

		if doy_arr_1[0]==doy_arr_2[0]:
#			del dataset['years']
#			del dataset['seq_months']
#			del dataset['seq_seasons']
			del dataset_2['years']
			del dataset_2['seq_months']
			del dataset_2['seq_seasons']

			print 'adjusting times in dataset_2 to those in dataset_1 to make diff possible.'

			dataset_2.coords['time']=(('time'),time_arr_1)

			dataset_2.coords['years']=(('time'),dataset.years.to_index())
			dataset_2.coords['seq_months']=(('time'),dataset.seq_months.to_index())
			dataset_2.coords['seq_seasons']=(('time'),dataset.seq_seasons.to_index())



		else:
			print 'FATAL - dates in DIFF file not lined up. First dayofyear in dataset ',doy_arr_1[0],' First dayofyear in dataset_2 ', doy_arr_2[0]
			sys.exit()

	dataset_diff=dataset-dataset_2

	return dataset_diff

def output_nc_file(dataset, field_name, model_params, output_dict):
    
    #create grid
    
    lons,lats,lonbs,latbs,nlon,nlat,nlonb,nlatb=cts.create_grid(output_dict['manual_grid_option'])
    
    if output_dict['is_thd']:
        p_full,p_half,npfull,nphalf=cts.create_pressures()
    else:
        p_full=None
        p_half=None
        npfull=None
        nphalf=None

    #create times
    try:
        output_dict['num_years']
    except KeyError:
        is_climatology=True
    else:
        if output_dict['num_years']==1:
            is_climatology=True
        else:
            is_climatology=False
        num_years=output_dict['num_years']
    
    time_arr,day_number,ntime,time_units=cts.create_time_arr(num_years,is_climatology,output_dict['time_spacing_days'])

    #Output it to a netcdf file. 
    file_name=output_dict['file_name']
    variable_name=output_dict['var_name']
    
    number_dict={}
    number_dict['nlat']=nlat
    number_dict['nlon']=nlon
    number_dict['nlatb']=nlatb
    number_dict['nlonb']=nlonb
    number_dict['npfull']=npfull
    number_dict['nphalf']=nphalf
    number_dict['ntime']=ntime
    
    data_out=dataset[field_name].load().data
    
    cts.output_to_file(data_out,lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,file_name,variable_name,number_dict)



