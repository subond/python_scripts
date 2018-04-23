import numpy as np
import filter as fl
import calendar_calc as cal
import cmip_time as cmt

def zonal_mean(array,zon_dim):

	zonal_mean=np.mean(array,axis=zon_dim)		

	return zonal_mean
	
def time_average(array,settings,time_dim, time_arr):
# 	time_average_settings = [time_averaged_plot,time_av_window_start,time_av_window_end]

	if not settings[0]:
		return array
			
	ndim=len(array.shape)
	start=settings[1]
	end=settings[2]
	averaging_type = settings[5]

	if averaging_type == 0:
		averaging_type = 'all_time'
	

	if averaging_type == 'all_time':
		if ndim==4 and time_dim==0:
			array_sub=array[start:end,:,:,:]
		elif ndim==3 and time_dim==0:
			array_sub=array[start:end,:,:]
		elif ndim==5 and time_dim==1:
			array_sub=array[:,start:end,:,:]
		else:
			print 'Something not right with time averaging'
			return -1

		time_mean=np.mean(array_sub,axis=time_dim)

	date_arr = cal.day_number_to_date(time_arr)
	subset_flag = False   #Initial assignment of subset_flag. 

	if averaging_type == 'djf':

		time_idx = (date_arr.month == 12) | (date_arr.month == 1) | (date_arr.month == 2)
		subset_flag = True

	if averaging_type == 'mam':

		time_idx = (date_arr.month == 3) | (date_arr.month == 4) | (date_arr.month == 5)
		subset_flag = True

	if averaging_type == 'jja':

		time_idx = (date_arr.month == 6) | (date_arr.month == 7) | (date_arr.month == 8)
		subset_flag = True

	if averaging_type == 'son':

		time_idx = (date_arr.month == 9) | (date_arr.month == 10) | (date_arr.month == 11)
		subset_flag = True

	if averaging_type == 'jan':

		time_idx = (date_arr.month == 1) 
		subset_flag = True

	if averaging_type == 'feb':

		time_idx = (date_arr.month == 2) 
		subset_flag = True

	if averaging_type == 'mar':

		time_idx = (date_arr.month == 3) 
		subset_flag = True

	if averaging_type == 'apr':

		time_idx = (date_arr.month == 4)
		subset_flag = True

	if averaging_type == 'may':

		time_idx = (date_arr.month == 5)
		subset_flag = True

	if averaging_type == 'jun':

		time_idx = (date_arr.month == 6)
		subset_flag = True

	if averaging_type == 'jul':

		time_idx = (date_arr.month == 7)
		subset_flag = True

	if averaging_type == 'aug':

		time_idx = (date_arr.month == 8)
		subset_flag = True

	if averaging_type == 'sep':

		time_idx = (date_arr.month == 9)
		subset_flag = True

	if averaging_type == 'oct':

		time_idx = (date_arr.month == 10)
		subset_flag = True

	if averaging_type == 'nov':

		time_idx = (date_arr.month == 11)
		subset_flag = True

	if averaging_type == 'dec':

		time_idx = (date_arr.month == 12)
		subset_flag = True

	if subset_flag:

		array_sub=array[time_idx,...]
		time_mean=np.mean(array_sub,axis=time_dim)

		

	return time_mean	


def perform_averaging(do_time_average, do_variance, do_time_filter, time_filter_bounds, time_average_settings, time_plot_number, units, input_data, time_arr):
	"Takes an input array with time and two spatial dimensions, and does the appropriate time averaging and / or filtering."

	if do_time_average:
		thd_data_temp=time_average(input_data,time_average_settings,0, time_arr)
		print 'done time average'
	elif np.logical_and(do_variance, do_time_filter):
                fs = fl.convert_frequencies_to_hertz( 1, 'days')
                lowcut = fl.convert_frequencies_to_hertz( time_filter_bounds[0], 'days')
                highcut = fl.convert_frequencies_to_hertz( time_filter_bounds[1], 'days')

                thd_data_temp = fl.butter_bandpass_filter(input_data, lowcut, highcut, fs, order=6, axis_to_filter=0)
                thd_data_temp=np.var(thd_data_temp, axis=0)
		print 'done variance with time filter'
	elif np.logical_and(do_variance, np.logical_not(do_time_filter)):

                thd_data_temp = input_data
                thd_data_temp=np.var(thd_data_temp, axis=0)

		print 'done variance without time filter'
	else: 
		thd_data_temp=np.squeeze(input_data[time_plot_number,...])
		print 'not done time average'

	thd_data=thd_data_temp
	del thd_data_temp

	if (units != 0.):
		thd_data = thd_data / units

	return thd_data




