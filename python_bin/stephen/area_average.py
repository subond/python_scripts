import set_and_get_params as sagp
import numpy as np

def area_average(dataset, variable_name, model_params, land_ocean_all='all', level=None, axis_in='time'):

	print 'performing area average on ',variable_name, 'of type ', land_ocean_all

	if (variable_name[0:9]=='hc_scaled'):
		data_to_average=dataset[variable_name[10:]]*dataset['ml_heat_cap']/model_params['delta_t']
	elif(variable_name[0:8]=='sigma_sb'):
		data_to_average=model_params['sigma_sb']*dataset[variable_name[9:]]**4.
	elif(level!=None):
		data_to_average=dataset[variable_name].sel(pfull=level, method='nearest')
	else:
		data_to_average=dataset[variable_name]

	try:
		grid_area=dataset['grid_cell_area']
	except KeyError:
		sagp.get_grid_sizes(dataset,model_params)
		grid_area=dataset['grid_cell_area']


	if(land_ocean_all == 'land'):
		scaled_grid_area=grid_area*(dataset['land'])

		multiplied=scaled_grid_area*data_to_average
		average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))

	elif(land_ocean_all == 'ocean'):
		scaled_grid_area=grid_area*(1.-dataset['land'])

		multiplied=scaled_grid_area*data_to_average
		average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))

	elif(land_ocean_all == 'all'):
		multiplied=grid_area*data_to_average
		average=multiplied.sum(('lat','lon'))/grid_area.sum(('lat','lon'))

	elif(land_ocean_all == 'qflux_area'):
		scaled_grid_area=grid_area*(dataset['qflux_area'])

		multiplied=scaled_grid_area*data_to_average
		average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))

	elif(land_ocean_all[3:] == 'eur'):
		scaled_grid_area=grid_area*(dataset[land_ocean_all])

		multiplied=scaled_grid_area*data_to_average
		average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))
	else:
		print 'invalid area-average option: ',land_ocean_all
		return

	new_var_name=variable_name+'_area_av_'+land_ocean_all

	dataset[new_var_name]=((axis_in), average)
	
def european_area_av(dataset, model_params, eur_area_av_input):

	variables_list=eur_area_av_input['variables_list']
	try:
        	levels_list  = eur_area_av_input['levels_list']
	except KeyError:
		levels_list  = None

	lats=dataset.lat
	lons=dataset.lon

	lon_array, lat_array = np.meshgrid(lons,lats)

	idx_nw_eur =     (45. <= lat_array) & (lat_array < 60.) & (-5. < lon_array) & (lon_array < 27.5)
	idx_nw_eur_neg = (45. <= lat_array) & (lat_array < 60.) & (np.mod(-5.,360.) < lon_array)

	idx_sw_eur =     (30. <= lat_array) & (lat_array < 45.) & (-5. < lon_array) & (lon_array < 27.5)
	idx_sw_eur_neg = (30. <= lat_array) & (lat_array < 45.) & (np.mod(-5.,360.) < lon_array)

	idx_ne_eur = (45. <= lat_array) & (lat_array < 60.) & (27.5 < lon_array) & (lon_array < 60.)
	idx_se_eur = (30. <= lat_array) & (lat_array < 45.) & (27.5 < lon_array) & (lon_array < 60.)


	land_nw_eur=np.zeros_like(dataset.land)
	land_sw_eur=np.zeros_like(dataset.land)

	land_ne_eur=np.zeros_like(dataset.land)
	land_se_eur=np.zeros_like(dataset.land)

	land_nw_eur[idx_nw_eur]=1.0
	land_nw_eur[idx_nw_eur_neg]=1.0

	land_sw_eur[idx_sw_eur]=1.0
	land_sw_eur[idx_sw_eur_neg]=1.0

	land_ne_eur[idx_ne_eur]=1.0
	land_se_eur[idx_se_eur]=1.0


	dataset['nw_eur']=(('lat','lon'), land_nw_eur)
	dataset['sw_eur']=(('lat','lon'), land_sw_eur)
	dataset['ne_eur']=(('lat','lon'), land_ne_eur)
	dataset['se_eur']=(('lat','lon'), land_se_eur)

	for i in range(np.shape(variables_list)[0]):
		var_name=variables_list[i]
		if levels_list!=None:
			level_in=levels_list[i]
		else:
			level_in=None

		area_average(dataset, var_name, model_params, land_ocean_all='nw_eur',level=level_in)
		area_average(dataset, var_name, model_params, land_ocean_all='sw_eur',level=level_in)
		area_average(dataset, var_name, model_params, land_ocean_all='ne_eur',level=level_in)
		area_average(dataset, var_name, model_params, land_ocean_all='se_eur',level=level_in)


def qflux_area_av(dataset, model_params, qflux_area_av_input):

	qflux_area=np.zeros_like(dataset.land)

	variables_list     = qflux_area_av_input['variables_list']

	warmpool_lat_centre= qflux_area_av_input['lat_centre']
	warmpool_lon_centre= qflux_area_av_input['lon_centre']

	warmpool_width     = qflux_area_av_input['width']
	warmpool_width_lon = qflux_area_av_input['width_lon']

	lats=dataset.lat
	lons=dataset.lon

	latbs=dataset.latb
	lonbs=dataset.lonb


	for j in np.arange(len(lats)):
	     lat = 0.5*(latbs[j+1] + latbs[j])
	     lat = (lat-warmpool_lat_centre)/warmpool_width
	     for i in np.arange(len(lons)):
		      lon = 0.5*(lonbs[i+1] + lonbs[i])
		      lon = (lon-warmpool_lon_centre)/warmpool_width_lon
		      if( lat**2.+lon**2. <= 1.0 ):
		          qflux_area[j,i]=1.0

	dataset['qflux_area']=(('lat','lon'), qflux_area)


	for i in range(np.shape(variables_list)[0]):
		var_name=variables_list[i]
		area_average(dataset, var_name, model_params, land_ocean_all='qflux_area')

import numpy as np
from cell_area import cell_area_all
import xarray as xarray

def planet_params_set(pref=1000.,kappa=2./7.,rm=6376.0e3, g=9.8, rct=6376.0e3, omega=7.292e-5):

	planet_params={}

	planet_params['pref']=pref
	planet_params['kappa']=kappa
	planet_params['rm']=rm #Planet radius in metres.
	planet_params['g']=g #Graviational acceleration
	planet_params['rct']=rct # Rgct is the R/g constant
	planet_params['rgct']=planet_params['rct']/planet_params['g']
	planet_params['omega'] = omega #Planet rotation rate omega

	return planet_params


def model_params_set(input_dir,res=42, delta_t=900., radius=6376.0e3):

	model_params={}

	model_params['input_dir']=input_dir
	model_params['res']=res
	model_params['delta_t']=delta_t
	model_params['planet_radius']=radius

	return model_params

def get_grid_sizes(dataset,model_params):

	area,x,y=cell_area_all(model_params['res'],model_params['input_dir'])

	area=(model_params['planet_radius']**2.)*area
	x=(model_params['planet_radius'])*x
	y=(model_params['planet_radius'])*y

	dataset['grid_cell_area']=(('lat','lon'),area)
	dataset['grid_cell_size_lon']=(('lat','lon'),x)
	dataset['grid_cell_size_lat']=(('lat','lon'),y)

def cfl_check(dataset,model_params):

	get_grid_sizes(dataset,model_params)
	dataset['cfl']=(('time','pfull','lat','lon'),np.abs(dataset['ucomp']*model_params['delta_t']*model_params['res']/model_params['planet_radius']))
	dataset['max_cfl']=(('time'),dataset['cfl'].max(dim=('pfull','lat','lon')))
#	dataset['lev_max_cfl']=(('time'),np.where(dataset.cfl == dataset.max_cfl)[1])

def run_analysis(analysis_list, dataset, model_params):

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
		


	return dataset



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

	
def rel_vort_calc(utheta,vtheta,lonar,latar,rm,n_arr,wholelat,wholelon):

# 	n_arr=[nlons,nlats,nlevs,ntheta,ntime_tot]

	nlons=n_arr[0]
	nlats=n_arr[1]
	nlevs=n_arr[2]	
	ntheta=n_arr[3]
	ntime_tot=n_arr[4]	

	# This bit does the absolute vorticity calculations.
	print 'Calculating Absolute Vorticity'
	dvdx=np.zeros((ntime_tot,ntheta,nlats,nlons))
	dudy=np.zeros((ntime_tot,ntheta,nlats,nlons))
	vort=np.zeros((ntime_tot,ntheta,nlats,nlons))

	pi=4.*np.arctan(1.)

	latarpi=((2*pi)/(360.0))*latar
	lonarpi=((2*pi)/(360.0))*lonar

# ; X-Gradient **********************************************************************************
# print,'running x grad'

	lats_diff=range(1,nlats-1)
	lons_diff=range(1,nlons-1)
	
	x_position=np.zeros((1,1,nlats,1))
	x_position[0,0,:,0]=rm*np.cos(range(nlats))
	x_position=x_position*np.ones((ntime_tot,ntheta,nlats,nlons))


	for nlon in lons_diff: 

		dvdx[:,:,:,nlon]=(vtheta[:,:,:,nlon+1]-vtheta[:,:,:,nlon-1])/((lonarpi[nlon+1]-lonarpi[nlon-1])*x_position[:,:,:,nlon])


	if wholelon:
		# x-gradient for first longitude (0)
		dvdx[:,:,:,0]=(vtheta[:,:,:,1]-vtheta[:,:,:,-1])/((lonarpi[2]-lonarpi[0])*x_position[:,:,:,0])

		# x-gradient for last longitude (n_elements(lonar)-1)
		dvdx[:,:,:,-1]=(vtheta[:,:,:,0]-vtheta[:,:,:,-2])/((lonarpi[2]-lonarpi[0])*x_position[:,:,:,-1])

	for nlat in lats_diff: 

		dudy[:,:,nlat,:]=(utheta[:,:,nlat+1,:]-utheta[:,:,nlat-1,:])/(rm*(latarpi[nlat+1]-latarpi[nlat-1]))-(utheta[:,:,nlat,:]*np.tan(latarpi[nlat])/rm)
		# For origin of second term see Middle Atmosphere Dynamics (Andrews et al) 1987 equation 3.8.4b. Minus sign to make sure that when we do dv/dx-du/dy we have the correct sign for the + u*tan(lat)/rm 

	vort=dvdx-dudy

	return vort
