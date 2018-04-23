"""
Evaluate and plot momentum budget at 150 hPa before and after onset

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams


def partition_advection(data, lons, lev=150):
    
    #First do uu terms
    uu_trans_dx = -86400. * gr.ddx( (data.ucomp_sq - data.ucomp**2).sel(pfull=lev) ) # <u'u'> = <uu> - <u><u>
    
    u = data.ucomp.sel(pfull=lev) # u
    u_dx = -86400. * gr.ddx( u )  # dudx
    
    u_dudx_zav = u.sel(lon=lons).mean('lon') * u_dx.sel(lon=lons).mean('lon') # [u][dudx]

    u_dudx_stat = u * u_dx - u_dudx_zav           # u*dudx* = [ududx] - [u][dudx]
    
    data['uu_trans_dx'] = (('lat','lon'), uu_trans_dx )	
    data['u_dudx_stat'] = (('lat','lon'), u_dudx_stat )	
    data['u_dudx_zav']  = (('lat'), u_dudx_zav )
    
    print 'uu terms done'
    
    #Next do uv terms
    uv_trans_dy = -86400. * gr.ddy( (data.ucomp_vcomp - data.ucomp * data.vcomp).sel(pfull=lev) , uv=True)

    v = data.vcomp.sel(pfull=lev).load() # v
    u_dy = -86400. * gr.ddy( u )  # dudy
    
    v_dudy_zav = v.sel(lon=lons).mean('lon') * u_dy.sel(lon=lons).mean('lon') # [v][dudy]
        
    v_dudy_stat = v * u_dy - v_dudy_zav           # v*dudy* = [vdudy] - [v][dudy]
        
    data['uv_trans_dy'] = (('lat','lon'), uv_trans_dy)	
    data['v_dudy_stat'] = (('lat','lon'), v_dudy_stat)	
    data['v_dudy_zav']  = (('lat'), v_dudy_zav )
    
    print 'uv terms done'
    
    #Finally do uw terms
    uw_trans_dp = -86400. * gr.ddp( (data.ucomp_omega - data.ucomp * data.omega))
    
    w = data.omega.sel(pfull=lev).load() # w
    u_dp = -86400. * (gr.ddp(data.ucomp)).sel(pfull=lev)  # dudp
    
    w_dudp_zav = w.sel(lon=lons).mean('lon') * u_dp.sel(lon=lons).mean('lon')
    w_dudp_stat = w * u_dp - w_dudp_zav
    
    data['uw_trans_dp'] = (('lat','lon'), uw_trans_dp.sel(pfull=lev))	
    data['w_dudp_stat'] = (('lat','lon'), w_dudp_stat)	
    data['w_dudp_zav']  = (('lat'), w_dudp_zav )	
    
    print 'uw terms done'
    
    
    
def mom_budg_latlon(run, t, lev=150, land=False,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 15, 10
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        
    if land:
        land_data = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')
    
    t_dates = data.xofyear[t:t+4]
    
    data_tmean = data.sel(xofyear=t_dates).mean('xofyear')
    
    #advective terms
    partition_advection(data_tmean, lons=lons, lev=lev)
    
    #Coriolis
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    fv = data_tmean.vcomp.sel(pfull=lev) * f * 86400.
    
    #Geopotential gradient
    dphidx = gr.ddx(data_tmean.height.sel(pfull=lev))
    dphidx = -86400. * 9.8 * dphidx

    mom_mean = data_tmean.u_dudx_zav + data_tmean.v_dudy_zav + data_tmean.w_dudp_zav
    mom_trans = data_tmean.uu_trans_dx + data_tmean.uv_trans_dy + data_tmean.uw_trans_dp
    mom_stat = data_tmean.u_dudx_stat + data_tmean.v_dudy_stat + data_tmean.w_dudp_stat
    
    mom_sum = fv + dphidx + mom_mean + mom_trans + mom_stat
    
    mom_mean = xr.DataArray( np.rollaxis(np.tile(mom_mean,[128,1]),1), mom_stat.coords)
    
    levels = np.arange(-30,30.1,5.)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    # Six subplots
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f1=fv.plot.contourf(ax=ax1, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)
    ax1.set_xlim(60,150)
    ax1.grid(True,linestyle=':')
    
    #Second plot
    mom_mean.plot.contourf(ax=ax2, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax2.grid(True,linestyle=':')
    ax2.set_ylim(-60,60)
    ax2.set_xlim(60,150)
    
    #Third plot
    dphidx.plot.contourf(ax=ax3, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax3.grid(True,linestyle=':')
    ax3.set_ylim(-60,60)
    ax3.set_xlim(60,150)
    
    #Fourth plot
    mom_trans.plot.contourf(ax=ax4, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax4.grid(True,linestyle=':')
    ax4.set_ylabel('Latitude')
    ax4.set_xlabel('Longitude')
    ax4.set_ylim(-60,60)
    ax4.set_xlim(60,150)
    ax4.set_xticks(range(60,151,15))
    
    #Fifth plot
    mom_stat.plot.contourf(ax=ax5, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax5.grid(True,linestyle=':')
    ax4.set_xlabel('Longitude')
    ax5.set_ylim(-60,60)
    ax5.set_xlim(60,150)
    ax5.set_xticks(range(60,151,15))
    
    #Sixth plot
    mom_sum.plot.contourf(ax=ax6, x='lon', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax6.grid(True,linestyle=':')
    ax4.set_xlabel('Longitude')
    ax6.set_ylim(-60,60)
    ax6.set_xlim(60,150)
    ax6.set_xticks(range(60,151,15))
    
    if land:
        for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
            land_data.land_mask.plot.contour(x='lon', y='lat', ax=ax, colors='w', levels=range(-1000,1002,1000), add_labels=False, add_colorbar=False)
    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
    #Colorbar
    cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
    cb1.set_label('$ms^{-1}day^{-1}$')
    
    figname = 'zon_mom_budg_ll_' + str(t) +'.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
        

mom_budg_latlon('ap_2', 30)
mom_budg_latlon('ap_2', 39)

mom_budg_latlon('full_qflux', 18, land=True, lonin=[60.,150.])
mom_budg_latlon('full_qflux', 39, land=True, lonin=[60.,150.])

mom_budg_latlon('flat_qflux', 18, land=True, lonin=[60.,150.])
mom_budg_latlon('flat_qflux', 44, land=True, lonin=[60.,150.])


