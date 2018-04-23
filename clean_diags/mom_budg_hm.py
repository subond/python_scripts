"""
Evaluate and plot momentum budget at 150 hPa

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

    u_dudx_stat = (u * u_dx).sel(lon=lons).mean('lon') - u_dudx_zav           # u*dudx* = [ududx] - [u][dudx]
    
    data['uu_trans_dx'] = (('xofyear','lat'), uu_trans_dx.sel(lon=lons).mean('lon') )	
    data['u_dudx_stat'] = (('xofyear','lat'), u_dudx_stat )	
    data['u_dudx_zav']  = (('xofyear','lat'), u_dudx_zav )
    
    print 'uu terms done'
    
    #Next do uv terms
    uv_trans_dy = -86400. * gr.ddy( (data.ucomp_vcomp - data.ucomp * data.vcomp).sel(pfull=lev) , uv=True)

    v = data.vcomp.sel(pfull=lev).load() # v
    u_dy = -86400. * gr.ddy( u )  # dudy
    
    v_dudy_zav = v.sel(lon=lons).mean('lon') * u_dy.sel(lon=lons).mean('lon') # [v][dudy]
        
    v_dudy_stat = (v * u_dy).sel(lon=lons).mean('lon') - v_dudy_zav           # v*dudy* = [vdudy] - [v][dudy]
        
    data['uv_trans_dy'] = (('xofyear','lat'), uv_trans_dy.sel(lon=lons).mean('lon'))	
    data['v_dudy_stat'] = (('xofyear','lat'), v_dudy_stat)	
    data['v_dudy_zav']  = (('xofyear','lat'), v_dudy_zav )
    
    print 'uv terms done'
    
    #Finally do uw terms
    uw_trans_dp = -86400. * gr.ddp( (data.ucomp_omega - data.ucomp * data.omega).sel(lon=lons).mean('lon') )
    
    w = data.omega.sel(pfull=lev).load() # w
    u_dp = -86400. * (gr.ddp(data.ucomp)).sel(pfull=lev)  # dudp
    
    w_dudp_zav = w.sel(lon=lons).mean('lon') * u_dp.sel(lon=lons).mean('lon')
    w_dudp_stat = (w * u_dp).sel(lon=lons).mean('lon') - w_dudp_zav
    
    data['uw_trans_dp'] = (('xofyear','lat'), uw_trans_dp.sel(pfull=lev))	
    data['w_dudp_stat'] = (('xofyear','lat'), w_dudp_stat)	
    data['w_dudp_zav']  = (('xofyear','lat'), w_dudp_zav )	
    
    print 'uw terms done'
    
    
    
def mom_budg_hm(run, months, lev=150, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 12, 10
    #rcParams['figure.figsize'] = 18, 10
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    #data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    #advective terms
    partition_advection(data, lons, lev=150)
    
    #Coriolis
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    fv = data.vcomp.sel(pfull=lev) * f * 86400.
    fv = fv.sel(lon=lons).mean('lon')
    
    #Geopotential gradient
    dphidx = gr.ddx(data.height.sel(pfull=lev))
    dphidx = -86400. * 9.8 * dphidx.sel(lon=lons).mean('lon')
        
    mom_mean = data.u_dudx_zav + data.v_dudy_zav + data.w_dudp_zav
    mom_trans = data.uu_trans_dx + data.uv_trans_dy + data.uw_trans_dp
    mom_stat = data.u_dudx_stat + data.v_dudy_stat + data.w_dudp_stat
    
    mom_sum = fv + dphidx + mom_mean + mom_trans + mom_stat
    
    levels = np.arange(-10,10.1,2.)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    # Six subplots
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f1=fv.plot.contourf(ax=ax1, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)
    ax1.grid(True,linestyle=':')
    
    #Second plot
    mom_mean.plot.contourf(ax=ax2, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax2.grid(True,linestyle=':')
    ax2.set_ylim(-60,60)
    
    #Third plot
    dphidx.plot.contourf(ax=ax3, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax3.grid(True,linestyle=':')
    ax3.set_ylim(-60,60)
    
    #Fourth plot
    mom_trans.plot.contourf(ax=ax4, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax4.grid(True,linestyle=':')
    ax4.set_ylabel('Latitude')
    ax4.set_xticks(tickspace)
    ax4.set_xticklabels(labels,rotation=25)
    ax4.set_ylim(-60,60)
    
    #Fifth plot
    mom_stat.plot.contourf(ax=ax5, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax5.grid(True,linestyle=':')
    ax5.set_xticks(tickspace)
    ax5.set_xticklabels(labels,rotation=25)
    ax5.set_ylim(-60,60)
    
    #Sixth plot
    mom_sum.plot.contourf(ax=ax6, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax6.grid(True,linestyle=':')
    ax6.set_xticks(tickspace)
    ax6.set_xticklabels(labels,rotation=25)
    ax6.set_ylim(-60,60)
    
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.3, wspace=0.1)
    #Colorbar
    cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
    #cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True,fraction=0.15, aspect=30)
    cb1.set_label('$ms^{-1}day^{-1}$')
    
    if lonin == [-1.,361.]:
        figname = 'zon_mom_budg.pdf'
    else:
        figname = 'zon_mom_budg_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
        
#mom_budg_hm('ap_20_qflux', [121,481])
#mom_budg_hm('ap_20_qflux', [121,481], lonin=[60.,150.])
#mom_budg_hm('am_qflux', [121,481])
#mom_budg_hm('am_qflux', [121,481], lonin=[60.,150.])
#mom_budg_hm('ap_2', [121,481])
#mom_budg_hm('full_qflux', [121,481])
mom_budg_hm('full_qflux', [121,481], lonin=[60.,90.])
mom_budg_hm('full_qflux', [121,481], lonin=[90.,120.])
mom_budg_hm('full_qflux', [121,481], lonin=[120.,150.])
#mom_budg_hm('flat_qflux', [121,481])
#mom_budg_hm('flat_qflux', [121,481], lonin=[60.,150.])


