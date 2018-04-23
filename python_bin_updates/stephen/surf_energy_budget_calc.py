import area_average as aav


def surf_energy_budget(dataset, model_params):

	aav.area_average(dataset, 't_surf', model_params)
	aav.area_average(dataset, 'sigma_sb_t_surf', model_params)
	aav.area_average(dataset, 'hc_scaled_delta_t_surf', model_params)
	aav.area_average(dataset, 'delta_t_surf', model_params)
	aav.area_average(dataset, 'flux_lhe', model_params)
	aav.area_average(dataset, 'flux_t', model_params)

	aav.area_average(dataset, 'flux_lw', model_params)
	aav.area_average(dataset, 'flux_sw', model_params)

	aav.area_average(dataset, 'sigma_sb_t_surf', model_params, land_ocean_all='land')
	aav.area_average(dataset, 't_surf', model_params, land_ocean_all='land')
	aav.area_average(dataset, 'hc_scaled_delta_t_surf', model_params, land_ocean_all='land')
	aav.area_average(dataset, 'delta_t_surf', model_params, land_ocean_all='land')
	aav.area_average(dataset, 'flux_lhe', model_params, land_ocean_all='land')
	aav.area_average(dataset, 'flux_t', model_params, land_ocean_all='land')

	aav.area_average(dataset, 'flux_lw', model_params, land_ocean_all='land')
	aav.area_average(dataset, 'flux_sw', model_params, land_ocean_all='land')

	rhs_land=dataset.flux_sw_area_av_land + dataset.flux_lw_area_av_land -dataset.sigma_sb_t_surf_area_av_land - dataset.flux_t_area_av_land - dataset.flux_lhe_area_av_land
	dataset['surf_energy_rhs_area_av_land']=(('time'),rhs_land)


	rhs_all=dataset.flux_sw_area_av_all + dataset.flux_lw_area_av_all -dataset.sigma_sb_t_surf_area_av_all - dataset.flux_t_area_av_all- dataset.flux_lhe_area_av_all
	dataset['surf_energy_rhs_area_av_all']=(('time'),rhs_all)


def convert_to_monthly(dataset, dataset_monthly, variable_name):

	data_to_monthly=dataset[variable_name].groupby('seq_months').mean('time')

	dataset_monthly[variable_name]=(('time','lat','lon'),data_to_monthly)