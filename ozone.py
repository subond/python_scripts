# -*- coding: utf-8 -*-s
import numpy as np
from stephen import create_timeseries as cts
import xarray as xr


def output_to_file(data,lats,lons,latbs,lonbs,p_full,p_half,file_name,variable_name,number_dict):

	output_file = Dataset(file_name, 'w', format='NETCDF3_CLASSIC')

	if p_full==None:
		is_thd=False
	else:
		is_thd=True


	lat = output_file.createDimension('lat', number_dict['nlat'])
	lon = output_file.createDimension('lon', number_dict['nlon'])

	latb = output_file.createDimension('latb', number_dict['nlatb'])
	lonb = output_file.createDimension('lonb', number_dict['nlonb'])

	if is_thd:
		pfull = output_file.createDimension('pfull', number_dict['npfull'])
		phalf = output_file.createDimension('phalf', number_dict['nphalf'])


	latitudes = output_file.createVariable('lat','d',('lat',))
	longitudes = output_file.createVariable('lon','d',('lon',))

	latitudebs = output_file.createVariable('latb','d',('latb',))
	longitudebs = output_file.createVariable('lonb','d',('lonb',))
	if is_thd:
		pfulls = output_file.createVariable('pfull','d',('pfull',))
		phalfs = output_file.createVariable('phalf','d',('phalf',))


	latitudes.units = 'degrees_N'.encode('utf-8')
	latitudes.cartesian_axis = 'Y'
	latitudes.edges = 'latb'
	latitudes.long_name = 'latitude'

	longitudes.units = 'degrees_E'.encode('utf-8')
	longitudes.cartesian_axis = 'X'
	longitudes.edges = 'lonb'
	longitudes.long_name = 'longitude'

	latitudebs.units = 'degrees_N'.encode('utf-8')
	latitudebs.cartesian_axis = 'Y'
	latitudebs.long_name = 'latitude edges'

	longitudebs.units = 'degrees_E'.encode('utf-8')
	longitudebs.cartesian_axis = 'X'
	longitudebs.long_name = 'longitude edges'

	if is_thd:
		pfulls.units = 'hPa'
		pfulls.cartesian_axis = 'Z'
		pfulls.positive = 'down'
		pfulls.long_name = 'full pressure level'

		phalfs.units = 'hPa'
		phalfs.cartesian_axis = 'Z'
		phalfs.positive = 'down'
		phalfs.long_name = 'half pressure level'

	if is_thd:
		output_array_netcdf = output_file.createVariable(variable_name,'f4',('pfull', 'lat','lon',))
	else:
		output_array_netcdf = output_file.createVariable(variable_name,'f4',('lat','lon',))

	latitudes[:] = lats
	longitudes[:] = lons

	latitudebs[:] = latbs
	longitudebs[:] = lonbs

	if is_thd:
		pfulls[:]     = p_full
		phalfs[:]     = p_half


	output_array_netcdf[:] = data

	output_file.close()







ozone_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/ozone_1990.nc'
ozone = xr.open_dataset( ozone_file, decode_times=False)


output_dict={'manual_grid_option':False, 'is_thd':True, 'num_years':1., 'time_spacing_days':12, 'file_name':'ozone_1990.nc', 'var_name':'ozone_1990'}

lons,lats,lonbs,latbs,nlon,nlat,nlonb,nlatb=cts.create_grid(output_dict['manual_grid_option'])

lons = ozone.lon.values

lonbs = ozone.lonb.values

nlon=lons.shape[0]

nlonb=lonbs.shape[0]
        
p_full=ozone.pfull.values
p_half=ozone.phalf.values
npfull=len(p_full)
nphalf=len(p_half)

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

data_out=ozone['ozone_1990'].load().data
print data_out.shape
data_out = data_out[0,:,:,:]

cts.output_to_file(data_out,lats,lons,latbs,lonbs,p_full,p_half,time_arr,time_units,file_name,variable_name,number_dict)

