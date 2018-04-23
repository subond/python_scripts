'''Function to write model variables to a netcdf file, e.g. for model input, 
or to allow the interpolator to be applied to new variables such as horizontal streamfunction.
Adapted from Stephens output functions'''

import numpy as np
import xarray as xr
import os
from netcdftime import utime
from netCDF4 import Dataset, date2num
from data_handling import FakeDT  # Mike's object for allowing time array to have additional attributes

GFDL_BASE = os.environ['GFDL_BASE']


def day_number_to_date(time_in, calendar_type = '360_day', units_in = 'days since 0001-01-01 00:00:00'):
    """
    Stephen: Aim is to make the time array have attributes like .month, or .year etc. This doesn't work with
    normal datetime objects, so Mike's FakeDT does this for you. First step is to turn input times
    into an array of datetime objects, and then FakeDT makes the array have the attributes of the
    elements themselves.
    """
    
    cdftime = utime(units_in, calendar = calendar_type)  # Create cdf time object
    time_in = cdftime.num2date(time_in)  # Use num2date method to convert time_in to this calendar
    cdftime = FakeDT( time_in, units = units_in,  
                 calendar = calendar_type)   #Create FakeDT object 
    return cdftime
    
    
    
def create_time_arr(is_climatology, num_per_year = 0., day_number = [0.]):
    
    if(is_climatology):
        num_days=360.
        day_number = np.linspace(0, num_days, num_per_year+1)[1:] - (num_days / (2.*num_per_year))
        time_units='days since 0000-01-01 00:00:00.0'
        print 'when creating a climatology file, the year of the time units must be zero. This is how the model knows it is a climatology.'
    else:
        # Use an input time array and set units as:
        time_units='days since 0001-01-01 00:00:00.0'
    print day_number
    time_arr = day_number_to_date(day_number, units_in = time_units)

    ntime=len(time_arr)
    
    return time_arr, time_units



def output_to_file(data, lats, lons, latbs, lonbs, p_full, p_half, time_arr, time_units, file_name, variable_name, number_dict):
    # Function to write data to a netcdf file 
    
    # Create output file to write to
    output_file = Dataset(file_name, 'w', format='NETCDF3_CLASSIC')
    
    # Identify whether we have pressure and time dimensions
    if p_full==None:
        pressure_dim = False
    else:
        pressure_dim = True
    
    if time_arr==None:
        time_dim = False
    else:
        time_dim = True
    
    # Create dimensions
    lat = output_file.createDimension('lat', number_dict['nlat'])
    lon = output_file.createDimension('lon', number_dict['nlon'])
    
    latb = output_file.createDimension('latb', number_dict['nlatb'])
    lonb = output_file.createDimension('lonb', number_dict['nlonb'])
    
    if pressure_dim:
        pfull = output_file.createDimension('pfull', number_dict['npfull'])
        phalf = output_file.createDimension('phalf', number_dict['nphalf'])
    
    if time_dim:
        time = output_file.createDimension('time', 0) #s Key point is to have the length of the time axis 0, or 'unlimited'. This seems necessary to get the code to run properly. 
    
    # Create variables for dimensions
    latitudes = output_file.createVariable('lat','d',('lat',))
    longitudes = output_file.createVariable('lon','d',('lon',))
    
    latitudebs = output_file.createVariable('latb','d',('latb',))
    longitudebs = output_file.createVariable('lonb','d',('lonb',))
    
    if pressure_dim:
        pfulls = output_file.createVariable('pfull','d',('pfull',))
        phalfs = output_file.createVariable('phalf','d',('phalf',))
    
    if time_dim:
        times = output_file.createVariable('time','d',('time',))
    
    #Create units for dimensions
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
    
    if pressure_dim:
        pfulls.units = 'hPa'
        pfulls.cartesian_axis = 'Z'
        pfulls.positive = 'down'
        pfulls.long_name = 'full pressure level'
        
        phalfs.units = 'hPa'
        phalfs.cartesian_axis = 'Z'
        phalfs.positive = 'down'
        phalfs.long_name = 'half pressure level'
    
    if time_dim:
        times.units = time_units
        times.calendar = 'THIRTY_DAY_MONTHS'
        times.calendar_type = 'THIRTY_DAY_MONTHS'
        times.cartesian_axis = 'T'
    
    # Create variable in output file
    if pressure_dim & time_dim:
        output_array_netcdf = output_file.createVariable(variable_name, 'f4', ('time', 'pfull', 'lat', 'lon',))
    elif pressure_dim:
        output_array_netcdf = output_file.createVariable(variable_name, 'f4', ('pfull', 'lat', 'lon',))
    elif time_dim:
        output_array_netcdf = output_file.createVariable(variable_name, 'f4', ('time', 'lat', 'lon',))
    else:
        output_array_netcdf = output_file.createVariable(variable_name, 'f4', ('lat', 'lon',))
    
    # Fill NetCDF with data
    latitudes[:] = lats
    longitudes[:] = lons
    
    latitudebs[:] = latbs
    longitudebs[:] = lonbs
    
    if pressure_dim:
        pfulls[:]     = p_full
        phalfs[:]     = p_half
    
    if time_dim:
        times[:]     = date2num(time_arr,units='days since 0001-01-01 00:00:00.0',calendar='360_day')
    
    print data.shape
    print output_array_netcdf
    
    output_array_netcdf[:] = data
    
    #Save and close
    output_file.close()






def xr_to_nc_file(dataset, field_name, output_dict):
    """Wrapper for output_to_file to write an xarray dataset to a netcdf file.
       Inputs:
       dataset: xarray dataset including the variable to plot
       field_name: name of field to plot
       output_dict: Dictionary detailing charactristics of output:
            is_thd: True if field to save is 3D
            notime_clim_tser: set to notime for no time axis, clim for climatology, tser for timeseries
            num_per_year: number of entries if climatology, used to calculate spacing
            file_name: name of file to save to
            var_name: variable name for output file
    """
    
    # Take grid from xarray dataset
    lons = dataset.lon.values
    lats = dataset.lat.values
    lonbs = dataset.lonb.values
    latbs = dataset.latb.values
    
    nlon=lons.shape[0]
    nlat=lats.shape[0]
    nlonb=lonbs.shape[0]
    nlatb=latbs.shape[0]
    
    
    # See if xarray dataset includes pressure:
    try:
        p_full = dataset.pfull.values
        p_half = dataset.phalf.values
    	npfull = len(p_full)
    	nphalf = len(p_half)
    except:
        p_full = None
        p_half = None
        npfull = None
        nphalf = None
    if not output_dict['is_thd']:
        p_full = None
        p_half = None
        npfull = None
        nphalf = None
        

    # Create times
    if output_dict['notime_clim_tser'] == 'notime':
        time_arr = None
        time_units = None
        
    elif output_dict['notime_clim_tser'] == 'clim':
        is_climatology = True
        try:
            output_dict['num_per_year']
        except:
            raise RuntimeError('Number of data points not provided for climatology')
        time_arr, time_units = create_time_arr( is_climatology, output_dict['num_per_year'] )
        
    elif output_dict['notime_clim_tser'] == 'tser':
        is_climatology = False
        day_number = dataset.time.values
        print day_number
        time_arr, time_units = create_time_arr( is_climatology, day_number = day_number )
        
    else:
        raise RuntimeError('Invalid input for notime_clim_tser')
        
    
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
    
    data_out=dataset[field_name].load().data
    
    output_to_file(data_out, lats, lons, latbs, lonbs, p_full, p_half, time_arr, time_units, file_name, variable_name, number_dict)

