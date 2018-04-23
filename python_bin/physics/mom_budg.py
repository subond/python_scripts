import numpy as np
from finite_difference import cfd
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import time_means

#Note: use interpolated data!

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
    

#def calc_mean_mom(u,v,w):
    #Calculate momentum advection by zonal and and time mean
    
#    a= 6376.0e3 #radius used in model
#    coslat = np.cos(u.lat * np.pi/180)
    
#    dudx = xr.DataArray( cfd(u.values, u.lon*np.pi/180,3),   [('xofyear', u.xofyear ), ('pfull', u.pfull ), ('lat', u.lat), ('lon', u.lon )])
#    ududx = u*dudx/coslat/a

#    dudy = xr.DataArray( cfd((u*coslat).values, u.lat*np.pi/180,2),   [('xofyear', u.xofyear ), ('pfull', u.pfull ), ('lat', u.lat), ('lon', u.lon )])
#    vdudy = v*dudy/coslat/a

#    dudp = xr.DataArray( cfd(u.values, u.pfull*100,1), [('xofyear', u.xofyear ), ('pfull', u.pfull ), ('lat', u.lat), ('lon', u.lon )])
#    wdudp = w*dudp
    
#    out = ududx + vdudy + wdudp
    
#    return out


def calc_mean_mom(uv,uw):
    #Calculate momentum advection by time and zonal mean
    
    a= 6376.0e3 #radius used in model
    coslat = np.cos(uv.lat * np.pi/180)
    
    uv = uv*coslat*coslat

    duvdy = xr.DataArray( cfd(uv.values ,uv.lat*np.pi/180,2),   [('xofyear', uv.xofyear ), ('pfull', uv.pfull ), ('lat', uv.lat)])
    duvdy = duvdy/coslat/coslat/a

    duwdp = xr.DataArray( cfd(uw.values,uv.pfull*100,1),   [('xofyear', uv.xofyear ), ('pfull', uv.pfull ), ('lat', uv.lat)])
    duwdp = duwdp
    
    out =  duvdy + duwdp
    
    return out
    


def calc_eddy_mom(uu,uv,uw):
    #Calculate momentum advection by eddies
    
    a= 6376.0e3 #radius used in model
    coslat = np.cos(uu.lat * np.pi/180)
    
    uv = uv*coslat*coslat
    
    duueddx = xr.DataArray( cfd(uu.values,uu.lon*np.pi/180,3),   [('xofyear', uu.xofyear ), ('pfull', uu.pfull ), ('lat', uu.lat), ('lon', uu.lon )])
    duueddx = duueddx/coslat/a

    duveddy = xr.DataArray( cfd(uv.values ,uu.lat*np.pi/180,2),   [('xofyear', uu.xofyear ), ('pfull', uu.pfull ), ('lat', uu.lat), ('lon', uu.lon )])
    duveddy = duveddy/coslat/coslat/a

    duweddp = xr.DataArray( cfd(uw.values,uu.pfull*100,1),   [('xofyear', uu.xofyear ), ('pfull', uu.pfull ), ('lat', uu.lat), ('lon', uu.lon )])
    duweddp = duweddp
    
    out = duueddx + duveddy + duweddp
    
    return out
    

def mombudg_fn(data):
    #Evaluate momentum budget
    
    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    f = 2 * omega * np.sin(data.lat *np.pi/180)
        
    uu_ms, uu_stat, uu_trans = eddy_decomp(data.ucomp, data.ucomp, data.ucomp_sq)
    uv_ms, uv_stat, uv_trans = eddy_decomp(data.ucomp, data.vcomp, data.ucomp_vcomp)
    uw_ms, uw_stat, uw_trans = eddy_decomp(data.ucomp, data.omega, data.ucomp_omega)
    print 'eddy decomposition done'
    
    mom_mean = calc_mean_mom( uv_ms, uw_ms)
    print 'mean advective terms done'
    mom_trans = calc_eddy_mom(uu_trans, uv_trans, uw_trans)
    print 'transient terms done'
    mom_stat = calc_eddy_mom(uu_stat, uv_stat, uw_stat)    
    print 'stationary terms done'
    
    dphidx = xr.DataArray( cfd(data.height.values,data.lon*np.pi/180,3),   [('xofyear', data.xofyear ), ('pfull', data.pfull ), ('lat', data.lat), ('lon', data.lon )])
    dphidx = -1*dphidx/coslat/a
    
    fv = data.vcomp*f
    
    
    data_out = xr.Dataset({'fv': fv, 'ddamp': data.dt_ug_diffusion, 'rdamp': data.udt_rdamp, 'dphidx': dphidx,
                     'mom_mean': mom_mean, 'mom_trans': mom_trans, 'mom_stat': mom_stat})

    return data_out




def mombudg_closure_fn(data):
    #Evaluate momentum budget
    
    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    mom_div = -1.*calc_eddy_mom(data.ucomp_sq, data.ucomp_vcomp, data.ucomp_omega)
    
    dphidx = xr.DataArray( cfd(data.height.values,data.lon*np.pi/180,3),   [('xofyear', data.xofyear ), ('pfull', data.pfull ), ('lat', data.lat), ('lon', data.lon )])
    dphidx = -1.*dphidx/coslat/a
    
    fv = data.vcomp*f
    
    mom_sum = (fv + dphidx*9.8 + mom_div)*10000.
    
    mom_sum[45,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-2.,2.1,0.2))
    plt.show()

    return 



    
    