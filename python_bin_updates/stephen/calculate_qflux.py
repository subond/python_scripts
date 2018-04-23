import numpy as np
from cell_area import cell_area_all
import xarray
from xarray import ufuncs as xruf
import time
import pdb
from scipy import stats
from mpl_toolkits.basemap import shiftgrid
import analyse as ana
import area_average as aav
import nc_file_io_xarray as io

def qflux_calc(dataset, model_params):
	"""A method for calculating seasonally-varying qfluxes, as described in Russell et al 1985 DOI:10.1016/0377-0265(85)90022-3"""
	upper_ocean_heat_content(dataset, model_params)
	net_surf_energy_flux(dataset, model_params)
	deep_ocean_heat_content(dataset, model_params)
	ocean_transport(dataset, model_params)

	output_dict={'manual_grid_option':False, 'is_thd':False, 'num_years':1., 'time_spacing_days':12, 'file_name':'seasonal_qflux_control_10m.nc', 'var_name':'ocean_qflux'}

	io.output_nc_file(dataset,'masked_ocean_transport', model_params, output_dict)
	

def upper_ocean_heat_content(dataset, model_params):
	"""Calculating upper-ocean heat content assuming a constant mixed layer depth, unlike Russel 1985, who have a seasonally-varying mixed layer depth. """
	print 'doing upper ocean heat content and rate of change calculation'
	sst_data=dataset['t_surf'].groupby('months').mean('time').load()
	weighted_sst_data=model_params['ocean_rho']*model_params['ocean_cp']*model_params['ml_depth']*sst_data*(1.0-dataset['land'])

	shape1=np.shape(weighted_sst_data)

	ny=shape1[1]
	nx=shape1[2]
	d_weighted_sst_data_dt=np.zeros_like(weighted_sst_data)
	delta_t=model_params['day_length']*30.

	for x in range(nx):
		for y in range(ny):
			d_weighted_sst_data_dt[:,y,x]=ana.time_gradient(weighted_sst_data[:,y,x], delta_t)

	dataset['d_weighted_sst_data_dt']=(('months_ax','lat','lon'), d_weighted_sst_data_dt)	


def net_surf_energy_flux(dataset, model_params):
	"""Calculates the net surface energy flux to be used in q-flux calcuation, but also calcualtes a scaling factor such that the annual average of the area-averaged surface flux is zero."""

	print 'doing net surf energy flux'

	aav.area_average(dataset, 'flux_sw', model_params, land_ocean_all='ocean')
	aav.area_average(dataset, 'flux_lw', model_params, land_ocean_all='ocean')
	aav.area_average(dataset, 'sigma_sb_t_surf', model_params, land_ocean_all='ocean')
	aav.area_average(dataset, 'flux_t', model_params, land_ocean_all='ocean')
	aav.area_average(dataset, 'flux_lhe', model_params, land_ocean_all='ocean')

	scaling_factor=(((dataset['sigma_sb_t_surf_area_av_ocean']+dataset['flux_t_area_av_ocean']+dataset['flux_lhe_area_av_ocean']-dataset['flux_lw_area_av_ocean'])/dataset['flux_sw_area_av_ocean']).groupby('months').mean('time')).mean('months')

	print 'using scale factor for SW of '+ str(scaling_factor)

	net_surf_energy_fl=(scaling_factor*dataset['flux_sw']+dataset['flux_lw']-(model_params['sigma_sb']*dataset['t_surf']**4.0)-dataset['flux_t']-dataset['flux_lhe'])*(1.0-dataset['land'])

	dataset['net_surf_energy_fl']=(('months_ax','lat','lon'), net_surf_energy_fl.groupby('months').mean('time'))
	aav.area_average(dataset, 'net_surf_energy_fl', model_params, land_ocean_all='ocean', axis_in='months_ax')
	
def deep_ocean_heat_content(dataset, model_params):

	print 'doing deep ocean heat content'
	aav.area_average(dataset, 'd_weighted_sst_data_dt', model_params, land_ocean_all='ocean', axis_in='months_ax')
	d_deep_ocean_dt=dataset['net_surf_energy_fl_area_av_ocean']-dataset['d_weighted_sst_data_dt_area_av_ocean']

	dataset['d_deep_ocean_dt']=(('months_ax'), d_deep_ocean_dt)	

def ocean_transport(dataset, model_params):

	ocean_transport=dataset['d_weighted_sst_data_dt']+dataset['d_deep_ocean_dt']-dataset['net_surf_energy_fl']
	masked_ocean_transport=ocean_transport*(1.0-dataset['land'])
	masked_land_transport=ocean_transport*(dataset['land'])


	dataset['ocean_transport']=(('months_ax','lat','lon'), ocean_transport)	
	dataset['masked_ocean_transport']=(('months_ax','lat','lon'), masked_ocean_transport)
	dataset['masked_land_transport']=(('months_ax','lat','lon'), masked_land_transport)	

if __name__ == "__main__":

	import nc_file_io_xarray as io
	import set_and_get_params as sagp

	input_dir='/scratch/rg419/GFDL_model/GFDLmoistModel/'
	base_dir='/scratch/rg419/Data_moist/'
	land_file='input/land.nc'
	base_exp_name='amip_10m/'
	exp_name='amip_10m/'

	start_file=121
	end_file=480
	land_present=True
	topo_present=False

	avg_or_daily='daily'

	model_params = sagp.model_params_set(input_dir, delta_t=720., ml_depth=10.)

	dataset, time_arr, size_list = io.read_data( base_dir,exp_name,start_file,end_file,avg_or_daily,topo_present)
    
	land_array, topo_array = io.read_land(input_dir,base_exp_name,land_present,True,size_list,land_file)
	dataset['land'] = (('lat','lon'),land_array)

	qflux_calc(dataset, model_params)


