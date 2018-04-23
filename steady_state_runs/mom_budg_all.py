# Make videos of the momentum budget for all terms and runs

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from physics import gradients as gr
from pylab import rcParams
import os 

rcParams['figure.figsize'] = 9, 12
rcParams['font.size'] = 18
rcParams['text.usetex'] = True


def mom_budg(run, lev=150):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
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
    
    data['uu_trans_dx'] = (('lat', 'lon'), uu_trans_dx )	
    data['u_dudx_cross1'] = (('lat', 'lon'), u_dudx_cross1 )	
    data['u_dudx_cross2'] = (('lat', 'lon'), u_dudx_cross2 )	
    data['u_dudx_stat'] = (('lat', 'lon'), u_dudx_stat )	
    data['u_dudx_zav']  = (('lat'), u_dudx_zav )
    
    
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
        
    data['uv_trans_dy'] = (('lat', 'lon'), uv_trans_dy)	
    data['v_dudy_cross1'] = (('lat', 'lon'), v_dudy_cross1 )	
    data['v_dudy_cross2'] = (('lat', 'lon'), v_dudy_cross2 )
    data['v_dudy_stat'] = (('lat', 'lon'), v_dudy_stat)	
    data['v_dudy_zav']  = (('lat'), v_dudy_zav )
    
    
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
    
    data['uw_trans_dp'] = (('lat', 'lon'), uw_trans_dp.sel(pfull=lev))	
    data['w_dudp_cross1'] = (('lat', 'lon'), w_dudp_cross1 )	
    data['w_dudp_cross2'] = (('lat', 'lon'), w_dudp_cross2 )
    data['w_dudp_stat'] = (('lat', 'lon'), w_dudp_stat)	
    data['w_dudp_zav']  = (('lat'), w_dudp_zav )	
    
    
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
    
    
    data['mom_mean'] = (('lat'), mom_mean)
    data['mom_cross'] = (('lat', 'lon'), mom_cross)
    data['mom_trans'] = (('lat', 'lon'), mom_trans)
    data['mom_stat'] = (('lat', 'lon'), mom_stat)
    data['fv_local'] = (('lat', 'lon'), fv_local)
    data['fv_ageo'] = (('lat', 'lon'), fv_ageo)
    data['fv_mean'] = (('lat'), fv_mean)
    data['dphidx'] = (('lat', 'lon'), dphidx)
    data['mom_sum'] = (('lat', 'lon'), mom_sum)
    
    return data


    
def subplot(data, ax_in, data_qflux):
    # Produce a subplot on a specified axis for a given pentad

    f1 = data.plot.contourf(ax=ax_in, x='lon', y='lat', extend='both', levels = np.arange(-40,41.1,4.), add_colorbar=False, add_labels=False)
    data_qflux.ocean_qflux.plot.contour(x='lon', y='lat', ax=ax_in, add_labels=False, colors='k', levels = np.arange(50.,350.,100.))
    ax_in.set_yticks(range(-60,61,30))
    ax_in.set_ylim(-60.,60.)
    ax_in.set_xticks(range(0,361,90))
    ax_in.grid(True,linestyle=':')
    
    return f1


def plot_var(var):
    # Plot up a given variable 
    
    plot_dir = '/scratch/rg419/plots/steady_state_runs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), 
     (ax10, ax11, ax12), (ax13, ax14, ax15), (ax16, ax17, ax18)) = plt.subplots(6, 3, sharex='col', sharey='row')
               
    plt.set_cmap('RdBu_r')
    
    ax_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17, ax18]
    lat_list = [0, 5, 10, 15, 20, 25]
    mag_list = [100, 200, 300]

    for i in range(18):
        mag = mag_list[np.mod(i,3)]
        lat = lat_list[i//3]
        #print mag, lat
        ax = ax_list[i]
    
        run = 'qflux_' + str(lat) + '_' + str(mag)
        data = mom_budg(run)
    
        qflux_name = os.environ["GFDL_BASE"] + 'exp/monsoon_ss_conts/input/qflux_95_lat_' + str(lat) + '_a_' + str(mag) + '.nc'
        data_qflux = xr.open_dataset(qflux_name)
        
        f1 = subplot(data[var], ax, data_qflux)


    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.05, hspace=0.25, wspace=0.15)

    cb1=f.colorbar(f1, ax=(ax_list), use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.07, aspect=30, shrink=0.5)

    figname = var + '_all.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    
plot_var('mom_cross')
plot_var('mom_trans')
plot_var('mom_stat')
plot_var('fv_local')
plot_var('dphidx')
plot_var('fv_ageo')