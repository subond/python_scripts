# Evaluate terms in the momentum budget for use by other modules (05/04/2018)

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from data_handling_updates import gradients as gr
from pylab import rcParams

rcParams['figure.figsize'] = 12, 4
rcParams['font.size'] = 16
rcParams['text.usetex'] = True



land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
land = xr.open_dataset( land_file)


def mom_budg_fun(run, lev=150):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/'+run+'.nc')
    
    #First do uu terms
    uu_trans_dx = -86400. * gr.ddx( (data.ucomp_sq - data.ucomp**2).sel(pfull=lev) ) # <u'u'> = <uu> - <u><u>
    
    u = data.ucomp.sel(pfull=lev) # u
    u_dx = -86400. * gr.ddx( u )  # dudx
    
    u_ed = u - u.mean('lon')
    u_dx_ed = u_dx - u_dx.mean('lon')
    
    u_dudx_zav = u.mean('lon') * u_dx.mean('lon') # [u][dudx], where brackets denote mean over all longitudes
    
    u_dudx_cross1 = u.mean('lon') * u_dx_ed # [u]dudx*
    
    u_dudx_cross2 = u_ed * u_dx.mean('lon') # u*[dudx]

    u_dudx_stat = u_ed * u_dx_ed         # u*dudx* 
    
    data['uu_trans_dx'] = (('xofyear','lat', 'lon'), uu_trans_dx )	
    data['u_dudx_cross1'] = (('xofyear','lat', 'lon'), u_dudx_cross1 )	
    data['u_dudx_cross2'] = (('xofyear','lat', 'lon'), u_dudx_cross2 )	
    data['u_dudx_stat'] = (('xofyear','lat', 'lon'), u_dudx_stat )	
    data['u_dudx_zav']  = (('xofyear','lat'), u_dudx_zav )
    
    
    #Next do uv terms
    uv_trans_dy = -86400. * gr.ddy( (data.ucomp_vcomp - data.ucomp * data.vcomp).sel(pfull=lev) , uv=True)

    v = data.vcomp.sel(pfull=lev).load() # v
    u_dy = -86400. * gr.ddy( u )  # dudy
    
    v_ed = v - v.mean('lon')
    u_dy_ed = u_dy - u_dy.mean('lon')
    
    v_dudy_zav = v.mean('lon') * u_dy.mean('lon') # [v][dudy]
    
    v_dudy_cross1 = v.mean('lon') * u_dy_ed # [v]dudy*
    
    v_dudy_cross2 = v_ed * u_dy.mean('lon')  # v*[dudy]
        
    v_dudy_stat = v_ed * u_dy_ed            # v*dudy* 
        
    data['uv_trans_dy'] = (('xofyear','lat', 'lon'), uv_trans_dy)	
    data['v_dudy_cross1'] = (('xofyear','lat', 'lon'), v_dudy_cross1 )	
    data['v_dudy_cross2'] = (('xofyear','lat', 'lon'), v_dudy_cross2 )
    data['v_dudy_stat'] = (('xofyear','lat', 'lon'), v_dudy_stat)	
    data['v_dudy_zav']  = (('xofyear','lat'), v_dudy_zav )
    
    
    #Finally do uw terms
    uw_trans_dp = -86400. * gr.ddp( (data.ucomp_omega - data.ucomp * data.omega))
    
    w = data.omega.sel(pfull=lev).load() # w
    u_dp = -86400. * (gr.ddp(data.ucomp)).sel(pfull=lev)  # dudp
    
    w_ed = w - w.mean('lon')
    u_dp_ed = u_dp - u_dp.mean('lon')
    
    w_dudp_zav = w.mean('lon') * u_dp.mean('lon') # [w][dudp]
    
    w_dudp_cross1 = w.mean('lon') * u_dp_ed # [w]dudp*
    
    w_dudp_cross2 = w_ed * u_dp.mean('lon') # w*[dudp]
    
    w_dudp_stat = w_ed * u_dp_ed         # w*dudp* 
    
    data['uw_trans_dp'] = (('xofyear','lat', 'lon'), uw_trans_dp.sel(pfull=lev))	
    data['w_dudp_cross1'] = (('xofyear','lat', 'lon'), w_dudp_cross1 )	
    data['w_dudp_cross2'] = (('xofyear','lat', 'lon'), w_dudp_cross2 )
    data['w_dudp_stat'] = (('xofyear','lat', 'lon'), w_dudp_stat)	
    data['w_dudp_zav']  = (('xofyear','lat'), w_dudp_zav )	
    
    
    #Coriolis
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    fv = data.vcomp.sel(pfull=lev) * f * 86400.
    fv_mean = fv.mean('lon')
    fv_local = fv - fv_mean
        
    #Geopotential gradient
    dphidx = gr.ddx(data.height.sel(pfull=lev))
    dphidx = -86400. * 9.8 * dphidx
    fv_ageo = fv_local + dphidx
        
    mom_mean = data.u_dudx_zav + data.v_dudy_zav + data.w_dudp_zav
    mom_cross = data.u_dudx_cross1 + data.v_dudy_cross1 + data.w_dudp_cross1 + data.u_dudx_cross2 + data.v_dudy_cross2 + data.w_dudp_cross2
    mom_trans = data.uu_trans_dx + data.uv_trans_dy + data.uw_trans_dp
    mom_stat = data.u_dudx_stat + data.v_dudy_stat + data.w_dudp_stat
    
    mom_sum = fv_local + fv_mean + dphidx + mom_mean + mom_trans + mom_stat + mom_cross
    
    
    data['mom_mean'] = (('xofyear','lat'), mom_mean)
    data['mom_cross'] = (('xofyear','lat', 'lon'), mom_cross)
    data['mom_trans'] = (('xofyear','lat', 'lon'), mom_trans)
    data['mom_stat'] = (('xofyear','lat', 'lon'), mom_stat)
    data['fv_local'] = (('xofyear','lat', 'lon'), fv_local)
    data['fv_ageo'] = (('xofyear','lat', 'lon'), fv_ageo)
    data['fv_mean'] = (('xofyear','lat'), fv_mean)
    data['dphidx'] = (('xofyear','lat', 'lon'), dphidx)
    data['mom_sum'] = (('xofyear','lat', 'lon'), mom_sum)
    
    return data
