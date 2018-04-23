import numpy as np
from cell_area import cell_area_all
import xarray
from xarray import ufuncs as xruf
import time
import pdb
import calculate_qflux as qfl
import set_and_get_params as sagp
import area_average as aav
import surf_energy_budget_calc as seb
import filter as fl

def run_analysis(analysis_list, dataset, dataset_monthly, model_params, eur_area_av_input=None, qflux_area_av_input=None):

	dataset_vars_read=['t_surf','temp']

	if analysis_list['ptemp']:
		print 'doing ptemp calc'
		dataset_vars_read, dataset=ana.pot_temp(dataset, lons, lats, pfull, time_arr, dataset_vars_read, planet_params) #Potential temperature calc. variable name 'ptemp'

	if analysis_list['brunt_vas']:
		print 'doing nsqd calc'
		nsqd_temp = ana.brunt_vas_freq(dataset, lons, lats, pfull, times, thd_vars_read , brunt_vas_level_1, brunt_vas_level_2, planet_params)
		dataset['nsqd']=nsqd_temp
		del nsqd_temp
		print 'doing eady calc'
		eady_temp = ana.eady_growth_rate(dataset, twd_data, lons, lats, pfull, times, thd_vars_read, twd_vars_read, brunt_vas_level_1, brunt_vas_level_2, planet_params)
		dataset['eady_gr']=eady_temp
		del eady_temp

	if analysis_list['pv']:
		print 'doing PV calc'
		#thetapre=np.array([50,100,150,200])
		thetapre=np.linspace(50,300,25)
		#Maybe mask this if less than 25?
		#print thetapre

		thd_vars_theta, dataset_theta, dataset, thd_vars_read, lons, lats =spv.pv_calc(dataset,thd_vars_read,thetapre,lons,lats,pfull,lons,lats, planet_params)

	if analysis_list['cfl']:
		print 'doing cfl calculation'
		cfl_check(dataset,model_params)

	if analysis_list['surf_energy_budget']:
		print 'doing surf_energy_budget'
		seb.surf_energy_budget(dataset, model_params)

	if analysis_list['european_area_av']:
		print 'doing european_area_av'
		aav.european_area_av(dataset, model_params, eur_area_av_input)

	if analysis_list['qflux_calc']:
		print 'doing qflux_calc'
		qfl.qflux_calc(dataset, model_params)

	if analysis_list['qflux_area_av']:
		print 'doing qflux_calc'
		aav.qflux_area_av(dataset, model_params, qflux_area_av_input)

	if analysis_list['merid_sf']:
		print 'doing meridsf_calc'
		merid_sf(dataset)

	return dataset

def cfl_check(dataset,model_params):

	sagp.get_grid_sizes(dataset,model_params)

	dataset_2=xruf.fabs(dataset['ucomp'])*model_params['delta_t']*model_params['res']/model_params['planet_radius']

	dataset['cfl']=(('time'),dataset_2.max(dim=('pfull','lat','lon')))

def time_gradient(data_in, delta_t):

	data_out=np.gradient(data_in, delta_t)

	return data_out	

def pot_temp(thd_data, lons, lats, pfull, times, vars_list, planet_params):

	pref=planet_params['pref']
	kappa=planet_params['kappa']

#	times_arr, pfull_arr, lat_arr, lon_arr = np.meshgrid(times,pfull,lats,lons,indexing='ij') #This creates a 4D array from pfull, so that the potential temperature calculations can be done without doing a loop.

#	del times_arr
#	del lat_arr
#	del lon_arr
	
	temp=thd_data['temp']
	
	pot_temp=temp*np.power((pref/pfull_arr),kappa)

#	pot_temp= pot_temp_individual( temp, pfull_arr, kappa, pref)

	vars_list.append('ptemp')
	thd_data['ptemp']=temp


	return vars_list, thd_data

def pot_temp_individual( temp, pfull, kappa = 2./7., pref=1000.):

	theta = temp * np.power((pref/pfull),kappa)

        return theta

	
def brunt_vas_freq(thd_data, lons, lats, pfull, times, vars_list, level_1, level_2, planet_params):

	grav = planet_params[1]

	temp_index=vars_list.index('temp')
	temp=np.squeeze(thd_data[temp_index,...])

	theta_index=vars_list.index('ptemp')
	theta=np.squeeze(thd_data[theta_index,...])

	height_index=vars_list.index('height')
	height=np.squeeze(thd_data[height_index,...])


	times_arr, pfull_arr, lat_arr, lon_arr = np.meshgrid(times,pfull,lats,lons,indexing='ij') #This creates a 4D array from pfull, so that the potential temperature calculations can be done without doing


#	p_centre = (pfull_arr[:,level_1,...] + pfull_arr[:,level_2,...])/2.
#	t_centre = (temp[:,level_1,...] + temp[:,level_2,...])/2.

        delta_z     = height[:,level_1,...] - height[:,level_2,...]
	delta_theta = theta [:,level_1,...] - theta [:,level_2,...] 

#	theta_env = pot_temp_individual( temp = t_centre , pfull = p_centre )
	theta_env = np.squeeze(theta [:,level_1,...] + theta [:,level_2,...] )/2.

	nsqd = np.squeeze((grav / theta_env) * (delta_theta / delta_z))


	return nsqd


def eady_growth_rate(thd_data, twd_data, lons, lats, pfull, times, vars_list, twd_vars_list, level_1, level_2, planet_params):


	times_arr, pfull_arr, lat_arr, lon_arr = np.meshgrid(times,pfull,lats,lons,indexing='ij')

	omega = planet_params[6]

	f_arr = 2. * omega * np.sin(np.pi*lat_arr[:,level_1,...]/180.)

	nsqd_index=twd_vars_list.index('nsqd')
	nsqd=np.squeeze(twd_data[nsqd_index,...])

	ucomp_index=vars_list.index('ucomp')
	ucomp=np.squeeze(thd_data[ucomp_index,...])

	height_index=vars_list.index('height')
	height=np.squeeze(thd_data[height_index,...])

        delta_z     = height[:,level_1,...] - height[:,level_2,...]
        delta_u     = ucomp[:,level_1,...]  - ucomp[:,level_2,...]


	eady_array = np.squeeze(0.31 * (np.abs(f_arr) / np.sqrt(nsqd) )* np.abs(delta_u / delta_z))


	unstable_nsqd_idx = (nsqd < 0.)
	eady_array[unstable_nsqd_idx] = 0. #if nsqd profile is unstable, i.e. negative, then we zero out the eady growth rate in this region. 
	

	return eady_array

def merid_sf(dataset, a=6371.0e3, g=9.8):
    """Calculate the mass streamfunction for the atmosphere.
    Based on a vertical integral of the meridional wind.
    Ref: Physics of Climate, Peixoto & Oort, 1992.  p158.
    `a` is the radius of the planet (default Earth 6371km).
    `g` is surface gravity (default Earth 9.8m/s^2).
    Returns an xarray DataArray of mass streamfunction.
    COPIED FROM https://github.com/ExeClim/ShareCode/blob/jpdev/execlim/analysis/mass_streamfunction.py on 16/05/16
    """
    vbar = dataset.vcomp.groupby('seasons').mean(('lon','time'))
    c = 2*np.pi*a*np.cos(vbar.lat*np.pi/180) / g
    # take a diff of half levels, and assign to pfull coordinates
    dp=xarray.DataArray(dataset.phalf.diff('phalf').values*100, coords=[('pfull', dataset.pfull)])
    merid_sf=c*np.cumsum(vbar*dp, axis=vbar.dims.index('pfull'))

    dataset['merid_sf']=(('seasons_ax','pfull','lat'),merid_st)




