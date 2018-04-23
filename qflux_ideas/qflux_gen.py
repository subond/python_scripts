#load in seasonally varying q fluxes
#average over given region for each month and normalise by integral so that at each longitude it sums to 0

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from stephen import aav, io, sagp

data = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/exp/q_flux_ap10/input/ocean_qflux.nc', decode_times=False)

lon_dic = {'all':range(0,128),
'africa':[i for i in range(len(data.lon)) if data.lon[i] >= 330. or data.lon[i] < 60.],
'asia':[i for i in range(len(data.lon)) if data.lon[i] >= 60. and data.lon[i] < 150.],
'pacific':[i for i in range(len(data.lon)) if data.lon[i] >= 150. and data.lon[i] < 240.],
'america':[i for i in range(len(data.lon)) if data.lon[i] >= 240. and data.lon[i] < 330.]}

def calc_qflux(qflux,lons,output_dict={0:0}):

    qflux['qflux_zav'] = (('time','lat','lon'), np.tile( np.mean(qflux.ocean_qflux.values[:,:,lon_dic[lons]],2,keepdims=True),[1,1,128]) )
    
    input_dir='/scratch/rg419/GFDL_model/GFDLmoistModel/'
    model_params = sagp.model_params_set(input_dir, delta_t=720., ml_depth=10.)
    aav.area_average(data, 'qflux_zav', model_params)
    
    qflux['qflux_out'] = (('time','lat','lon'), qflux.qflux_zav - qflux.qflux_zav_area_av_all)
    
    #io.output_nc_file(qflux,'qflux_out', model_params, output_dict)
    
    return qflux

def integrate_lat(qflux,lons):
    a= 6376.0e3 #radius used in model
    plot_dir='/scratch/rg419/plots/qflux_run_analysis/heat_trans/'

    qflux = calc_qflux(data,lons)
    
    qflux.coords['dlat'] = ('lat',np.diff(qflux.latb))
    qflux.coords['dlon'] = ('lon',np.diff(qflux.lonb))
    integrand = qflux.qflux_out * qflux.dlat*np.pi/180. * a * a * np.cos(qflux.lat * np.pi/180.) * qflux.dlon*np.pi/180.
    
    integrand_lon = (integrand).sum(('lon'))

    qflux['heat_trans'] = (('time','lat'), np.cumsum(integrand_lon.values, axis=1) )
    
    qflux.heat_trans.mean('time').plot()
    qflux.heat_trans[[0,1,11],:].mean('time').plot()
    qflux.heat_trans[2:4,:].mean('time').plot()
    qflux.heat_trans[5:7,:].mean('time').plot()
    qflux.heat_trans[8:10,:].mean('time').plot()
    plt.legend(['Year mean','DJF','MAM','JJA','SON'],loc=2)
    plt.ylim(-2.5e16,2.5e16)
    plt.savefig(plot_dir+'ocean_heat_trans_'+lons+'.png')
    plt.clf()
    
    qflux.heat_trans.plot(x='time',y='lat',vmin=-2.5e16,vmax=2.5e16)
    plt.savefig(plot_dir+'oht_hm_'+lons+'.png')
    plt.clf()

#integrate_lat(data,'asia')
#integrate_lat(data,'africa')
#integrate_lat(data,'pacific')
#integrate_lat(data,'america')
integrate_lat(data,'all')

output_dict={'manual_grid_option':False, 'is_thd':False, 'num_years':1., 'time_spacing_days':12, 'file_name':'ocean_qflux_all.nc', 'var_name':'ocean_qflux'}
#qflux_all = calc_qflux(data,'all',output_dict)

output_dict={'manual_grid_option':False, 'is_thd':False, 'num_years':1., 'time_spacing_days':12, 'file_name':'ocean_qflux_asia.nc', 'var_name':'ocean_qflux'}
#qflux_asia = calc_qflux(data,asia,output_dict)

output_dict={'manual_grid_option':False, 'is_thd':False, 'num_years':1., 'time_spacing_days':12, 'file_name':'ocean_qflux_africa.nc', 'var_name':'ocean_qflux'}
#qflux_africa = calc_qflux(data,africa,output_dict)

output_dict={'manual_grid_option':False, 'is_thd':False, 'num_years':1., 'time_spacing_days':12, 'file_name':'ocean_qflux_pacific.nc', 'var_name':'ocean_qflux'}
#qflux_pacific = calc_qflux(data,pacific,output_dict)

output_dict={'manual_grid_option':False, 'is_thd':False, 'num_years':1., 'time_spacing_days':12, 'file_name':'ocean_qflux_america.nc', 'var_name':'ocean_qflux'}
#qflux_america = calc_qflux(data,america,output_dict)



