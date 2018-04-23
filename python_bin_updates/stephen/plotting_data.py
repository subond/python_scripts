import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
import averaging as av
import os
import filter as fl

def run_plotting(dataset, plot_list):

	for plot_var in plot_list.iterkeys():

		pt_set=dict(plot_list[plot_var]) #Making pt_set a new copy of the dictionary each time, as otherwise it retains previous entried, which is a classic python of array assignment only being a label for one piece of memory, rather than permenant assignment. 
		try:
			plot_data = dataset[pt_set['var_name']].pfull
			is_thd=True
		except AttributeError:
			is_thd=False

		if is_thd:
			print 'plotting 3D data for '+pt_set['var_name']+' '+pt_set['plot_type']+' '+str(plot_var)
			plot_thd(dataset, pt_set)
		else:
			print 'plotting 2D data for '+pt_set['var_name']+' '+pt_set['plot_type']+' '+str(plot_var)
			plot_twd(dataset, pt_set)


def plot_thd(thd_data, pt_set):

	pt_set['exp_name']=thd_data.exp_name

	if(pt_set['plot_type']=='latlon'):
		try:
			plot_data = thd_data[pt_set['var_name']].sel(pfull=pt_set['plev'], method='nearest').groupby(pt_set['group_type']).mean(pt_set['av_type'])
		except AttributeError:
			plot_data = thd_data[pt_set['var_name']].sel(pfull=pt_set['plev'], method='nearest').mean(pt_set['av_type'])

		ntime=get_ntime(plot_data,pt_set)
		pt_set['ntime']=ntime

		land_array = thd_data['land']
		pt_set=plotting_presets(plot_data,pt_set)
		plot_latlon(plot_data,pt_set, land_array)


	if(pt_set['plot_type']=='latp'):
		try:
			plot_data = thd_data[pt_set['var_name']].groupby(pt_set['group_type']).mean(pt_set['av_type'])
		except AttributeError:
			plot_data = thd_data[pt_set['var_name']].mean(pt_set['av_type'])

		ntime=get_ntime(plot_data,pt_set)
		pt_set['ntime']=ntime

		nsub_rows,nsub_cols=get_row_col(ntime)

		plot_data.plot.contourf(x='lat', y='pfull', col=pt_set['group_type'], col_wrap=nsub_cols, levels=25)
		plt.ylim(thd_data.pfull.max(), thd_data.pfull.min())
		save_figure(pt_set)


def plot_twd(twd_data, pt_set):

	pt_set['exp_name']=twd_data.exp_name

	try:
		plot_data = twd_data[pt_set['var_name']].groupby(pt_set['group_type']).mean(pt_set['av_type'])
	except AttributeError:
		plot_data = twd_data[pt_set['var_name']].mean(pt_set['av_type'])
	land_array = twd_data['land']
	pt_set=plotting_presets(plot_data,pt_set)

	ntime=get_ntime(plot_data,pt_set)
	pt_set['ntime']=ntime


	if(pt_set['plot_type']=='latlon'):
		plot_latlon(plot_data,pt_set, land_array)

	if(pt_set['plot_type']=='zmlat'):
#		plot_data.plot.line(col=pt_set['group_type'], col_wrap=4, levels=25)
		plt.plot(plot_data[0,:],plot_data.lat)
		save_figure(pt_set)



def plot_latlon(plot_data, pt_set, land_array):

	ntime=pt_set['ntime']

	lon_centre=0.
	lat_centre=0.
	m = Basemap(projection='kav7',lat_0=lat_centre,lon_0=lon_centre)
	plot_data_temp,lons = shiftgrid(180.,plot_data,plot_data.lon,start=False)
	if np.max(land_array) > 0.:
		land_array_temp,lons = shiftgrid(180.,land_array,land_array.lon,start=False)
		plot_land = True
	else:
		plot_land = False

	lats = plot_data.lat.to_index()
	lon, lat = np.meshgrid(lons, lats)
	fig = plt.figure()
	levels=create_levels(plot_data, pt_set )



	for time_dim in np.arange(ntime):

		nsub_rows,nsub_cols=get_row_col(ntime)

		plt.subplot(nsub_rows,nsub_cols,time_dim+1)
		xi, yi = m(lon, lat)
		cs = m.contourf(xi,yi,plot_data_temp[time_dim,...],levels, cmap=plt.get_cmap(pt_set['cmap_name']))
	        m.drawparallels(np.arange(-90.,120.,30.))
        	m.drawmeridians(np.arange(0.,360.,60.))
		plt.title(pt_set['group_type']+'=%d' % time_dim)
		if plot_land:
			land = m.contour(xi,yi,land_array_temp)

	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	fig.colorbar(cs, cax=cbar_ax,extend='both')

	if(pt_set['exp_name'][0:5]=='diff_'):
		exp_title='diff'
	else:
		exp_title=pt_set['exp_name']
	title=pt_set['var_name'] +' in '+exp_title+' averaged over '+pt_set['av_type'] + ' when grouped by '+pt_set['group_type']
	title = title.translate(None, '!@#$/')

	fig.suptitle(title)
	save_figure(pt_set)

def create_levels(plot_data, pt_set):

	centre_on_zero=pt_set['centre_on_zero']
	presc_min_max=pt_set['presc_min_max']
	nlevels=pt_set['nlevels']
	goes_below_zero=pt_set['goes_below_zero']
	
	
	if goes_below_zero and nlevels%2==0: #This is so that IF you have data that goes below zero, and you prescribe an EVEN number of levels, this will make sure zero is included.
			nlevels=nlevels+1
		
	if presc_min_max:
		presc_min=presc_min_max[0]
		presc_max=presc_min_max[1]	
	elif goes_below_zero and centre_on_zero: 
		plot_data_max_in=np.max(plot_data)
		plot_data_min_in=np.min(plot_data)
		
		plot_data_max, plot_data_min = round_max_min(plot_data_max_in, plot_data_min_in)

		presc_min=-np.max([np.absolute(plot_data_min),np.absolute(plot_data_max)])
		presc_max= np.max([np.absolute(plot_data_min),np.absolute(plot_data_max)])
	else:
		presc_max_in=np.max(plot_data)
		presc_min_in=np.min(plot_data)
		presc_max, presc_min = round_max_min(presc_max_in, presc_min_in)


	
	levels = np.linspace(presc_min,presc_max, nlevels)
	return levels



def plotting_presets(plot_data,pt_set_in):

	pt_set_out=pt_set_in


	try:
		pt_set_out['presc_min_max']
	except KeyError:
		pt_set_out['presc_min_max'] = False

	try:
		pt_set_out['nlevels']
	except KeyError:
		pt_set_out['nlevels'] = 25

	try:
		pt_set_out['goes_below_zero']
	except KeyError:
		try:
			pt_set_out['goes_below_zero'] = plot_data.valid_range[0] < 0.
		except AttributeError:
			pt_set_out['goes_below_zero'] = pt_set_in['exp_name'][0:4] == 'diff' #one reason it wouldn't have valid_range defined is if it's a difference, meaning it must go below zero. (Could later have arrays that I create that don't have a valid range...)

	try:
		pt_set_out['centre_on_zero']
	except KeyError:
		pt_set_out['centre_on_zero'] = pt_set_out['goes_below_zero']

	try:
		pt_set_out['cmap_name']
	except KeyError:
		if(pt_set_out['goes_below_zero']):
			pt_set_out['cmap_name'] = 'RdBu_r'
		else:
			pt_set_out['cmap_name'] = 'RdBu_r'

	return pt_set_out

def round_max_min(max_in, min_in):
	if(max_in >= 1.0 or min_in <= -1.0):
		max_out = np.ceil(max_in)
		min_out = np.floor(min_in)
	else:
		max_out=max_in
		min_out=min_in
	return max_out, min_out


def get_row_col(ntotal):
	""" get row col finds the optimal arrangement of the number of plots available. iterate)factor(ntotal) works in all cases except ntotal=1 and ntotal = 2. We therefore run a while loop here that doesn't run iterate factor when these conditions are met. """
	factor_1=1
	factor_2=ntotal
	nitt=ntotal

	while(factor_1==1 and factor_2>2):
		factor_1,factor_2=iterate_factor(nitt)
		nitt=ntotal+1

	return factor_1, factor_2

def iterate_factor(ntotal):
	""" iterate_facor returns the pair of factors of ntotal that are closest together. I.e. iterate_factor(12)=3,4 (not 1,12 or 2,6).   """
	val_c=int(np.ceil(np.sqrt(ntotal)))
	val_f=int(np.floor(np.sqrt(ntotal)))
	if(val_c == val_f):
		factor_1=val_c
		factor_2=val_c
	else:
		fac=val_f
		mod_res=np.mod(ntotal,fac)
		while(mod_res !=0):
			fac=fac-1
			mod_res=np.mod(ntotal,fac)

		factor_1=int(fac)
		factor_2=int(ntotal/fac)

	return factor_1, factor_2	

def get_ntime(plot_data,pt_set):

	if(pt_set['group_type']=='seasons'):
		ntime=np.shape(plot_data.seasons.to_index())[0]
	elif(pt_set['group_type']=='years'):
		ntime=np.shape(plot_data.years.to_index())[0]
	elif(pt_set['group_type']=='months'):
		ntime=np.shape(plot_data.months.to_index())[0]
	elif(pt_set['group_type']=='all_time'):
		ntime=1

	return ntime

def save_figure(pt_set):

	if len(pt_set['av_type'][0])==1:
		av_name=pt_set['av_type'] 
	else:
		av_name='_'.join(pt_set['av_type']) 

	name = pt_set['var_name']+'_'+pt_set['exp_name']+'_'+av_name+'_'+str(pt_set['ntime'])+'_'+pt_set['group_type']
	name = name.translate(None, '!@#$/')

	exp_name=pt_set['exp_name']

	if(exp_name[0:5]=='diff_'):
		directory='/scratch/sit204/plots/diffs/'+exp_name[5:].translate(None, '!@#$/')+'/'
	else:
		directory='/scratch/sit204/plots/exps/'+exp_name.translate(None, '!@#$/')+'/'
	if not os.path.exists(directory):
		os.makedirs(directory)



	plt.savefig(directory+name+'.pdf', bbox_inches='tight')


