#Create ozone single time file.
import xarray as xr
import numpy as np
from data_handling import time_means
from finite_difference import cfd
import matplotlib.pyplot as plt

data = time_means('aquaplanet_10m_test', [121,145], filename='plev_pentad', timeav='month')

def eddy_decomp(a,b,ab):
    #Decompose fields into mean, stationary and transient eddies
    
    a_zmean = a.mean(('lon'))
    a_zed = a - a_zmean
    b_zmean = b.mean(('lon'))
    b_zed = b - b_zmean
    
    ab_ms = a_zmean*b_zmean
    ab_stat = a_zmean*b_zed + a_zed*b_zmean + a_zed*b_zed
    ab_trans = ab - ab_stat - ab_ms
    
    return ab_ms, ab_stat, ab_trans
    


def calc_mean_heat(vt,wt):
    #Calculate momentum advection by time and zonal mean
    
    a= 6376.0e3 #radius used in model
    coslat = np.cos(vt.lat * np.pi/180)
    
    vt = vt*coslat

    dvtdy = xr.DataArray( cfd(vt.values ,vt.lat*np.pi/180,2),   [('xofyear', vt.xofyear ), ('pfull', vt.pfull ), ('lat', vt.lat)])
    dvtdy = duvdy/coslat/a

    dwtdp = xr.DataArray( cfd(uw.values,vt.pfull*100,1),   [('xofyear', vt.xofyear ), ('pfull', vt.pfull ), ('lat', vt.lat)])
    dwtdp = dwtdp
    
    out =  dvtdy + dwtdp
    
    return out
    


def calc_eddy_heat(ut,vt,wt):
    #Calculate momentum advection by eddies
    
    a= 6376.0e3 #radius used in model
    coslat = np.cos(ut.lat * np.pi/180.)
    
    vt = vt*coslat
    
    duteddx = xr.DataArray( cfd(ut.values,vt.lon*np.pi/180.,3),   [('xofyear', vt.xofyear ), ('pfull', vt.pfull ), ('lat', vt.lat), ('lon', vt.lon )])
    duteddx = duteddx/coslat/a

    dvteddy = xr.DataArray( cfd(vt.values ,vt.lat*np.pi/180.,2),   [('xofyear', vt.xofyear ), ('pfull', vt.pfull ), ('lat', vt.lat), ('lon', vt.lon )])
    dvteddy = dvteddy/coslat/a

    dwteddp = xr.DataArray( cfd(wt.values,vt.pfull*100.,1),   [('xofyear', vt.xofyear ), ('pfull', vt.pfull ), ('lat', vt.lat), ('lon', vt.lon )])
    dwteddp = dwteddp
    
    out = duteddx + dvteddy + dwteddp
    
    return dvteddy, dwteddp
    

def heatbudg_fn(data):
    #Evaluate momentum budget
    
    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    f = 2 * omega * np.sin(data.lat *np.pi/180)
        
    ut_ms, ut_stat, ut_trans = eddy_decomp(data.ucomp, data.temp, data.ucomp_temp)
    vt_ms, vt_stat, vt_trans = eddy_decomp(data.vcomp, data.temp, data.vcomp_temp)
    wt_ms, wt_stat, wt_trans = eddy_decomp(data.omega, data.temp, data.omega_temp)
    print 'eddy decomposition done'
    
    heat_mean = calc_mean_heat( vt_ms, wt_ms)
    print 'mean advective terms done'
    heat_trans = calc_eddy_heat(ut_trans, vt_trans, wt_trans)
    print 'transient terms done'
    heat_stat = calc_eddy_heat(ut_stat, vt_stat, wt_stat)    
    print 'stationary terms done'
    
    pressure_term = data.omega * data.temp /data.pfull
        
    
    data_out = xr.Dataset({'cnvdamp': data.dt_tg_convection, 'difdamp': data.dt_tg_diffusion, 'cnddamp': data.dt_tg_condensation, 'raddamp': data.tdt_rad,
                     'heat_mean': mom_mean, 'heat_trans': mom_trans, 'heat_stat': mom_stat})

    return data_out




def heatbudg_closure_fn(data):
    #Evaluate momentum budget
    
    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    heat_div_y, heat_div_p = calc_eddy_heat(data.ucomp_temp, data.vcomp_temp, data.omega_temp)
    heat_div_y = heat_div_y*-86400.
    heat_div_p = heat_div_p*-86400.
    heat_diab = (data.dt_tg_convection + data.dt_tg_condensation + data.dt_tg_diffusion + data.tdt_rad)*86400.
    pressure_term = 2./7. * data.omega * data.temp /data.pfull*86400./100.
    
    heat_sum = (heat_diab +  pressure_term + heat_div_y + heat_div_p)
    
    #plt.figure(0)
    #heat_div[5,2,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-10,11,1))
    #plt.figure(0)
    #(data.dt_tg_convection[5,:,:,:]*86400.).mean('lon').plot.contourf(x='lat', y='pfull',yincrease=False, levels=np.arange(-6,6,0.5))
    #plt.figure(1)
    #(data.dt_tg_condensation[5,:,:,:]*86400.).mean('lon').plot.contourf(x='lat', y='pfull', yincrease=False, levels=np.arange(-6,6,0.5))
    #plt.figure(2)
    #(data.dt_tg_diffusion[5,:,:,:]*86400.).mean('lon').plot.contourf(x='lat', y='pfull', yincrease=False, levels=np.arange(-6,6,0.5))
    #plt.figure(3)
    #(data.tdt_rad[5,:,:,:]*86400.).mean('lon').plot.contourf(x='lat', y='pfull', yincrease=False, levels=np.arange(-6,6,0.5))
    #plt.figure(3)
    #heat_diab.mean(('xofyear','lon')).plot.contourf(x='lat', y='pfull', yincrease=False, levels=np.arange(-6,6,0.5))
    #plt.figure(4)
    #pressure_term.mean(('xofyear','lon')).plot.contourf(x='lat', y='pfull', yincrease=False, levels=np.arange(-6,6,0.5))
    plt.figure(5)
    heat_div_y.mean(('xofyear','lon')).plot.contourf(x='lat', y='pfull', yincrease=False, levels=np.arange(-6,6,0.5))
    plt.figure(6)
    heat_div_p.mean(('xofyear','lon')).plot.contourf(x='lat', y='pfull', yincrease=False, levels=np.arange(-6,6,0.5))
    
    
    #plt.figure(2)
    #heat_sum[5,2,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-10,11,1))
    #plt.figure(1)
    #heat_sum[5,2:17,:,:].mean('lon').plot.contourf(x='lat', y='pfull')
    plt.show()

    return 

heatbudg_closure_fn(data)