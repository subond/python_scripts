import numpy as np
import xarray as xr
from stephen import ana, aav, io, sagp


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
    sst_data=dataset.t_surf.load()
    weighted_sst_data=model_params['ocean_rho']*model_params['ocean_cp']*model_params['ml_depth']*sst_data*(1.0-dataset.land)
    
    shape1=np.shape(weighted_sst_data)
    
    ny=shape1[1]
    nx=shape1[2]
    d_weighted_sst_data_dt=np.zeros_like(weighted_sst_data)
    delta_t=model_params['day_length']*30.
    
    for x in range(nx):
        for y in range(ny):
            d_weighted_sst_data_dt[:,y,x]=ana.time_gradient(weighted_sst_data[:,y,x], delta_t)
    
    dataset['d_weighted_sst_data_dt']=(('xofyear','lat','lon'), d_weighted_sst_data_dt)


def net_surf_energy_flux(dataset, model_params):
    """Calculates the net surface energy flux to be used in q-flux calcuation, but also calcualtes a scaling factor such that the annual average of the area-averaged surface flux is zero."""
    
    print 'doing net surf energy flux'
    aav.area_average(dataset, 'flux_sw', model_params, land_ocean_all='ocean', axis_in='xofyear')
    aav.area_average(dataset, 'flux_lw', model_params, land_ocean_all='ocean', axis_in='xofyear')
    aav.area_average(dataset, 'sigma_sb_t_surf', model_params, land_ocean_all='ocean', axis_in='xofyear')
    aav.area_average(dataset, 'flux_t', model_params, land_ocean_all='ocean', axis_in='xofyear')
    aav.area_average(dataset, 'flux_lhe', model_params, land_ocean_all='ocean', axis_in='xofyear')
    
    scaling_factor=((dataset['sigma_sb_t_surf_area_av_ocean']+dataset['flux_t_area_av_ocean']+dataset['flux_lhe_area_av_ocean']-dataset['flux_lw_area_av_ocean'])/dataset['flux_sw_area_av_ocean']).mean('xofyear')
    
    print 'using scale factor for SW of '+ str(scaling_factor)
    
    net_surf_energy_fl=(scaling_factor*dataset['flux_sw']+dataset['flux_lw']-(model_params['sigma_sb']*dataset['t_surf']**4.0)-dataset['flux_t']-dataset['flux_lhe'])*(1.0-dataset['land'])
    
    dataset['net_surf_energy_fl']=(('xofyear','lat','lon'), net_surf_energy_fl)
    aav.area_average(dataset, 'net_surf_energy_fl', model_params, land_ocean_all='ocean', axis_in='xofyear')
    
	
def deep_ocean_heat_content(dataset, model_params):

	print 'doing deep ocean heat content'
	aav.area_average(dataset, 'd_weighted_sst_data_dt', model_params, land_ocean_all='ocean', axis_in='xofyear')
	d_deep_ocean_dt=dataset['net_surf_energy_fl_area_av_ocean']-dataset['d_weighted_sst_data_dt_area_av_ocean']

	dataset['d_deep_ocean_dt']=(('xofyear'), d_deep_ocean_dt)	

def ocean_transport(dataset, model_params):

	ocean_transport=dataset['d_weighted_sst_data_dt']+dataset['d_deep_ocean_dt']-dataset['net_surf_energy_fl']
	masked_ocean_transport=ocean_transport*(1.0-dataset['land'])
	masked_land_transport=ocean_transport*(dataset['land'])


	dataset['ocean_transport']=(('xofyear','lat','lon'), ocean_transport)	
	dataset['masked_ocean_transport']=(('xofyear','lat','lon'), masked_ocean_transport)
	dataset['masked_land_transport']=(('xofyear','lat','lon'), masked_land_transport)	

if __name__ == "__main__":
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/control_fixed_sst_2mtest.nc')
    data.coords['month'] = (data.xofyear-1)//6 +1
    data = data.groupby('month').mean(('xofyear'))   
        
    land = xr.open_dataset('/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    #land = xr.open_dataset('/scratch/rg419/python_scripts/land_era/land_era_no_america.nc')
    data['land'] = (('lat','lon'),land.land_mask.values)
    model_params = sagp.model_params_set('/scratch/rg419/Isca/', delta_t=720., ml_depth=10.)
    qflux_calc(data, model_params)


