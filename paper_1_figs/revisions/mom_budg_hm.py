"""
Evaluate and plot momentum budget at 150 hPa - redo in line with reviewer 1s suggestions, see how it looks

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling_updates import time_means, month_dic, gradients as gr
import sh
#from physics import gradients as gr
from pylab import rcParams


def partition_advection(data, lons, lev=150):
    
    #First do uu terms
    uu_trans_dx = -86400. * gr.ddx( (data.ucomp_sq - data.ucomp**2).sel(pfull=lev) ) # <u'u'> = <uu> - <u><u>
    
    u = data.ucomp.sel(pfull=lev) # u
    u_dx = -86400. * gr.ddx( u )  # dudx
    
    u_ed = u - u.mean('lon')
    u_dx_ed = u_dx - u_dx.mean('lon')
    
    u_dudx_zav = u.mean('lon') * u_dx.mean('lon') # [u][dudx], where brackets denote mean over all longitudes
    
    u_dudx_cross1 = (u.mean('lon') * u_dx_ed).sel(lon=lons).mean('lon') # [u]dudx*
    
    u_dudx_cross2 = (u_ed * u_dx.mean('lon')).sel(lon=lons).mean('lon') # u*[dudx]

    u_dudx_stat = (u_ed * u_dx_ed).sel(lon=lons).mean('lon')          # u*dudx* 
    
    data['uu_trans_dx'] = (('xofyear','lat'), uu_trans_dx.sel(lon=lons).mean('lon') )	
    data['u_dudx_cross1'] = (('xofyear','lat'), u_dudx_cross1 )	
    data['u_dudx_cross2'] = (('xofyear','lat'), u_dudx_cross2 )	
    data['u_dudx_stat'] = (('xofyear','lat'), u_dudx_stat )	
    data['u_dudx_zav']  = (('xofyear','lat'), u_dudx_zav )
    
    print('uu terms done')
    
    #Next do uv terms
    uv_trans_dy = -86400. * gr.ddy( (data.ucomp_vcomp - data.ucomp * data.vcomp).sel(pfull=lev) , uv=True)

    v = data.vcomp.sel(pfull=lev).load() # v
    u_dy = -86400. * gr.ddy( u )  # dudy
    
    v_ed = v - v.mean('lon')
    u_dy_ed = u_dy - u_dy.mean('lon')
    
    v_dudy_zav = v.mean('lon') * u_dy.mean('lon') # [v][dudy]
    
    v_dudy_cross1 = (v.mean('lon') * u_dy_ed).sel(lon=lons).mean('lon') # [v]dudy*
    
    v_dudy_cross2 = (v_ed * u_dy.mean('lon')).sel(lon=lons).mean('lon') # v*[dudy]
        
    v_dudy_stat = (v_ed * u_dy_ed).sel(lon=lons).mean('lon')           # v*dudy* 
        
    data['uv_trans_dy'] = (('xofyear','lat'), uv_trans_dy.sel(lon=lons).mean('lon'))	
    data['v_dudy_cross1'] = (('xofyear','lat'), v_dudy_cross1 )	
    data['v_dudy_cross2'] = (('xofyear','lat'), v_dudy_cross2 )
    data['v_dudy_stat'] = (('xofyear','lat'), v_dudy_stat)	
    data['v_dudy_zav']  = (('xofyear','lat'), v_dudy_zav )
    
    print('uv terms done')
    
    #Finally do uw terms
    uw_trans_dp = -86400. * gr.ddp( (data.ucomp_omega - data.ucomp * data.omega).sel(lon=lons).mean('lon') )
    
    w = data.omega.sel(pfull=lev).load() # w
    u_dp = -86400. * (gr.ddp(data.ucomp)).sel(pfull=lev)  # dudp
    
    w_ed = w - w.mean('lon')
    u_dp_ed = u_dp - u_dp.mean('lon')
    
    w_dudp_zav = w.mean('lon') * u_dp.mean('lon') # [w][dudp]
    
    w_dudp_cross1 = (w.mean('lon') * u_dp_ed).sel(lon=lons).mean('lon') # [w]dudp*
    
    w_dudp_cross2 = (w_ed * u_dp.mean('lon')).sel(lon=lons).mean('lon') # w*[dudp]
    
    w_dudp_stat = (w_ed * u_dp_ed).sel(lon=lons).mean('lon')         # w*dudp* 
    
    data['uw_trans_dp'] = (('xofyear','lat'), uw_trans_dp.sel(pfull=lev))	
    data['w_dudp_cross1'] = (('xofyear','lat'), w_dudp_cross1 )	
    data['w_dudp_cross2'] = (('xofyear','lat'), w_dudp_cross2 )
    data['w_dudp_stat'] = (('xofyear','lat'), w_dudp_stat)	
    data['w_dudp_zav']  = (('xofyear','lat'), w_dudp_zav )	
    
    print('uw terms done')
    
    
    
def mom_budg_hm(run, lev=150, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 12, 8
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/paper_1_figs/revisions/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/'+run+'.nc')
    
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
    fv_mean = fv.mean('lon')
    fv_local = (fv - fv_mean).sel(lon=lons).mean('lon')
    
    try:
        totp = ((data.precipitation)*86400.).sel(lon=lons).mean('lon')
    except:
        totp = ((data.convection_rain + data.condensation_rain)*86400.).sel(lon=lons).mean('lon')
    #abs_vort = (data.vor + f).sel(lon=lons).mean('lon')*86400.
    
    #Geopotential gradient
    dphidx = gr.ddx(data.height.sel(pfull=lev))
    dphidx = -86400. * 9.8 * dphidx.sel(lon=lons).mean('lon')
    fv_ageo = fv_local + dphidx
        
    mom_mean = data.u_dudx_zav + data.v_dudy_zav + data.w_dudp_zav
    mom_cross = data.u_dudx_cross1 + data.v_dudy_cross1 + data.w_dudp_cross1 + data.u_dudx_cross2 + data.v_dudy_cross2 + data.w_dudp_cross2
    mom_cross1 = data.u_dudx_cross1 + data.v_dudy_cross1 + data.w_dudp_cross1 
    mom_cross2 = data.u_dudx_cross2 + data.v_dudy_cross2 + data.w_dudp_cross2
    mom_trans = data.uu_trans_dx + data.uv_trans_dy + data.uw_trans_dp
    mom_stat = data.u_dudx_stat + data.v_dudy_stat + data.w_dudp_stat
    
    mom_crossu  = data.u_dudx_cross1 + data.u_dudx_cross2
    mom_crossv  = data.v_dudy_cross1 + data.v_dudy_cross2
    mom_crossw  = data.w_dudp_cross1 + data.w_dudp_cross2
    
    #print mom_stat
    
    mom_sum = fv_local + fv_mean + dphidx + mom_mean + mom_trans + mom_stat + mom_cross
    
    levels = np.arange(-20,21.1,2.)
    
    mn_dic = month_dic(1)
    tickspace = list(range(13,72,18))
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    # Six subplots
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f1=fv_mean.plot.contourf(ax=ax1, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    totp.plot.contour(ax=ax1, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax1.set_ylabel('Latitude')
    ax1.set_title('Zonal mean Coriolis', fontsize=17)
    ax1.set_ylim(-60,60)
    ax1.grid(True,linestyle=':')
    ax1.set_yticks(np.arange(-60.,61.,30.))
    ax1.text(-15, 60, 'a)')
    
    #Second plot
    mom_mean.plot.contourf(ax=ax2, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    #ax2.contour(data.xofyear, data.lat, abs_vort.sel(pfull=lev).T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2)
    #abs_vort.sel(pfull=lev).plot.contour(ax=ax2, x='xofyear', y='lat', extend='both', levels=np.arange(-2.,6.,2.), add_colorbar=False, add_labels=False, colors='k', linewidths=2)
    totp.plot.contour(ax=ax2, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax2.grid(True,linestyle=':')
    ax2.set_title('Mean state advection', fontsize=17)
    ax2.set_ylim(-60,60)
    ax2.text(-5, 60, 'b)')
    
    #Third plot
    mom_sum.plot.contourf(ax=ax3, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    totp.plot.contour(ax=ax3, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax3.grid(True,linestyle=':')
    ax3.set_ylim(-60,60)
    ax3.set_title('Residual', fontsize=17)
    ax3.text(-5, 60, 'c)')
    
    #Fourth plot
    mom_trans.plot.contourf(ax=ax4, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    totp.plot.contour(ax=ax4, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax4.grid(True,linestyle=':')
    ax4.set_ylabel('Latitude')
    ax4.set_ylim(-60,60)
    ax4.set_title('Transient eddy flux conv.', fontsize=17)
    ax4.set_yticks(np.arange(-60.,61.,30.))
    ax4.text(-15, 60, 'd)')
    
    #Fifth plot
    mom_stat.plot.contourf(ax=ax5, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    totp.plot.contour(ax=ax5, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax5.grid(True,linestyle=':')
    ax5.set_title('Stat. eddy flux conv.', fontsize=17)
    ax5.set_ylim(-60,60)
    ax5.text(-5, 60, 'e)')
    
    #Sixth plot
    mom_cross.plot.contourf(ax=ax6, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    totp.plot.contour(ax=ax6, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax6.grid(True,linestyle=':')
    ax6.set_ylim(-60,60)
    ax6.set_title('Stat. eddy cross terms', fontsize=17)
    ax6.text(-5, 60, 'f)')
    
    #Fourth plot
    dphidx.plot.contourf(ax=ax7, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    totp.plot.contour(ax=ax7, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax7.grid(True,linestyle=':')
    ax7.set_ylabel('Latitude')
    ax7.set_xticks(tickspace)
    ax7.set_xticklabels(labels,rotation=25)
    ax7.set_ylim(-60,60)
    ax7.set_title('Geopotential gradient', fontsize=17)
    ax7.set_yticks(np.arange(-60.,61.,30.))
    ax7.text(-15, 60, 'g)')
    
    #Fifth plot
    fv_local.plot.contourf(ax=ax8, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    totp.plot.contour(ax=ax8, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax8.grid(True,linestyle=':')
    ax8.set_xticks(tickspace)
    ax8.set_xticklabels(labels,rotation=25)
    ax8.set_title('Stat. eddy Coriolis', fontsize=17)
    ax8.set_ylim(-60,60)
    ax8.text(-5, 60, 'h)')
    
    #Sixth plot
    fv_ageo.plot.contourf(ax=ax9, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    totp.plot.contour(ax=ax9, x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
    ax9.grid(True,linestyle=':')
    ax9.set_xticks(tickspace)
    ax9.set_xticklabels(labels,rotation=25)
    ax9.set_ylim(-60,60)
    ax9.set_title('Ageostrophic stat. eddy Coriolis', fontsize=17)
    ax9.text(-5, 60, 'i)')
    
    
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.25, wspace=0.12)
    #Colorbar
    cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    #cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True,fraction=0.15, aspect=30)
    #cb1.set_label('$ms^{-1}day^{-1}$')
    
    if lonin == [-1.,361.]:
        figname = 'zon_mom_budg_' +run+ '_rev.pdf'
    else:
        figname = 'zon_mom_budg_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '_rev.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()

#mom_budg_hm('zs_sst')
#mom_budg_hm('ap10_qflux')
#mom_budg_hm('ap10_co2')

#mom_budg_hm('ap_2')
#mom_budg_hm('full_qflux')
mom_budg_hm('full_qflux', lonin=[60.,150.])
#mom_budg_hm('full_qflux', lonin=[345.,45.])
mom_budg_hm('flat_qflux', lonin=[60.,150.])
mom_budg_hm('am_qflux', lonin=[60.,150.])
mom_budg_hm('ap_20_qflux', lonin=[60.,150.])


