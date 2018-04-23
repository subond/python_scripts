import numpy as np
import pdb
import cell_area as carea

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


def model_params_set(input_dir,res=42, delta_t=900., radius=6376.0e3, day_length=86400., ocean_rho=1.035e3, ocean_cp=3989.24495292815, ml_depth=10., sigma_sb=5.6734e-8):

	model_params={}

	model_params['input_dir']=input_dir
	model_params['res']=res
	model_params['delta_t']=delta_t
	model_params['planet_radius']=radius
	model_params['sigma_sb']=sigma_sb
	model_params['day_length']=day_length
	model_params['ocean_rho']=ocean_rho
	model_params['ocean_cp']=ocean_cp
	model_params['ml_depth']=ml_depth



	return model_params

def get_grid_sizes(dataset,model_params):

	area,x,y=carea.cell_area_all(model_params['res'],model_params['input_dir'])

	area=(model_params['planet_radius']**2.)*area
	x=(model_params['planet_radius'])*x
	y=(model_params['planet_radius'])*y

	dataset['grid_cell_area']=(('lat','lon'),area)
	dataset['grid_cell_size_lon']=(('lat','lon'),x)
	dataset['grid_cell_size_lat']=(('lat','lon'),y)

