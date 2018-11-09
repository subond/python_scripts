"""
Evaluate and plot meridional momentum budget at 150 hPa 

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling_updates import month_dic, gradients as gr
import sh
from pylab import rcParams


def partition_advection(data, lons, lev=150):
    
    #First do uu terms
    uv_trans_dx = -86400. * gr.ddx( (data.ucomp_vcomp - data.ucomp*data.vcomp).sel(pfull=lev) ) # <u'u'> = <uu> - <u><u>
    
    u = data.ucomp.sel(pfull=lev) # u
    v = data.vcomp.sel(pfull=lev) # v
    
    v_dx = -86400. * gr.ddx( v )  # dvdx
    
    u_ed = u - u.mean('lon')
    v_dx_ed = v_dx - v_dx.mean('lon')
    
    u_dvdx_zav = u.mean('lon') * v_dx.mean('lon') # [u][dudx], where brackets denote mean over all longitudes
    
    u_dvdx_cross1 = (u.mean('lon') * v_dx_ed).sel(lon=lons).mean('lon') # [u]dudx*
    
    u_dvdx_cross2 = (u_ed * v_dx.mean('lon')).sel(lon=lons).mean('lon') # u*[dudx]

    u_dvdx_stat = (u_ed * v_dx_ed).sel(lon=lons).mean('lon')          # u*dudx* 
    
    data['uv_trans_dx'] = (('xofyear','lat'), uv_trans_dx.sel(lon=lons).mean('lon') )	
    data['u_dvdx_cross1'] = (('xofyear','lat'), u_dvdx_cross1 )	
    data['u_dvdx_cross2'] = (('xofyear','lat'), u_dvdx_cross2 )	
    data['u_dvdx_stat'] = (('xofyear','lat'), u_dvdx_stat )	
    data['u_dvdx_zav']  = (('xofyear','lat'), u_dvdx_zav )
    
    print('uv terms done')
    
    #Next do vv terms
    vv_trans_dy = -86400. * gr.ddy( (data.vcomp_sq - data.vcomp * data.vcomp).sel(pfull=lev) , uv=True)

    v_dy = -86400. * gr.ddy( v)  # dvdy
    
    v_ed = v - v.mean('lon')
    v_dy_ed = v_dy - v_dy.mean('lon')
    
    v_dvdy_zav = v.mean('lon') * v_dy.mean('lon') # [v][dudy]
    
    v_dvdy_cross1 = (v.mean('lon') * v_dy_ed).sel(lon=lons).mean('lon') # [v]dvdy*
    
    v_dvdy_cross2 = (v_ed * v_dy.mean('lon')).sel(lon=lons).mean('lon') # v*[dvdy]
        
    v_dvdy_stat = (v_ed * v_dy_ed).sel(lon=lons).mean('lon')           # v*dvdy* 
        
    data['vv_trans_dy'] = (('xofyear','lat'), vv_trans_dy.sel(lon=lons).mean('lon'))	
    data['v_dvdy_cross1'] = (('xofyear','lat'), v_dvdy_cross1 )	
    data['v_dvdy_cross2'] = (('xofyear','lat'), v_dvdy_cross2 )
    data['v_dvdy_stat'] = (('xofyear','lat'), v_dvdy_stat)	
    data['v_dvdy_zav']  = (('xofyear','lat'), v_dvdy_zav )
    
    print('vv terms done')
    
    #Finally do vw terms
    vw_trans_dp = -86400. * gr.ddp( (data.vcomp_omega - data.vcomp * data.omega).sel(lon=lons).mean('lon') )
    
    w = data.omega.sel(pfull=lev) # w
    v_dp = -86400. * (gr.ddp(data.vcomp)).sel(pfull=lev)  # dvdp
    
    w_ed = w - w.mean('lon')
    v_dp_ed = v_dp - v_dp.mean('lon')
    
    w_dvdp_zav = w.mean('lon') * v_dp.mean('lon') # [w][dvdp]
    
    w_dvdp_cross1 = (w.mean('lon') * v_dp_ed).sel(lon=lons).mean('lon') # [w]dvdp*
    
    w_dvdp_cross2 = (w_ed * v_dp.mean('lon')).sel(lon=lons).mean('lon') # w*[dvdp]
    
    w_dvdp_stat = (w_ed * v_dp_ed).sel(lon=lons).mean('lon')         # w*dvdp* 
    
    data['vw_trans_dp'] = (('xofyear','lat'), vw_trans_dp.sel(pfull=lev))	
    data['w_dvdp_cross1'] = (('xofyear','lat'), w_dvdp_cross1 )	
    data['w_dvdp_cross2'] = (('xofyear','lat'), w_dvdp_cross2 )
    data['w_dvdp_stat'] = (('xofyear','lat'), w_dvdp_stat)	
    data['w_dvdp_zav']  = (('xofyear','lat'), w_dvdp_zav )	
    
    print('vw terms done')
    
    
    
def mom_budg_hm(run, lev=150, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.], plot_precip=True, rot_fac=1.):
    
    rcParams['figure.figsize'] = 12, 6
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/mom_budg_merid/'
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
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    fu = -1. * data.ucomp.sel(pfull=lev) * f * 86400.
    fu_mean = fu.mean('lon')
    fu_local = (fu - fu_mean).sel(lon=lons).mean('lon')
    
    if plot_precip:
        try:
            totp = ((data.precipitation)*86400.).sel(lon=lons).mean('lon')
        except:
            totp = ((data.convection_rain + data.condensation_rain)*86400.).sel(lon=lons).mean('lon')
    #abs_vort = (data.vor + f).sel(lon=lons).mean('lon')*86400.
    
    #Geopotential gradient
    dphidy = gr.ddy(data.height.sel(pfull=lev), vector=False)
    dphidy = -86400. * 9.8 * dphidy.sel(lon=lons).mean('lon')
    fu_ageo = fu_local + fu_mean + dphidy
    
    dvdt = gr.ddt(data.vcomp.sel(pfull=lev)).sel(lon=lons).mean('lon')*86400.
    
    metric = (-86400.* data.ucomp_sq * np.tan(data.lat * np.pi/180.)/6376.0e3).sel(pfull=lev).sel(lon=lons).mean('lon')
    vert_term = (-86400.*data.vcomp_omega/6376.0e3).sel(pfull=lev).sel(lon=lons).mean('lon')
    
    mom_mean = data.u_dvdx_zav + data.v_dvdy_zav + data.w_dvdp_zav
    mom_cross = data.u_dvdx_cross1 + data.v_dvdy_cross1 + data.w_dvdp_cross1 + data.u_dvdx_cross2 + data.v_dvdy_cross2 + data.w_dvdp_cross2
    mom_cross1 = data.u_dvdx_cross1 + data.v_dvdy_cross1 + data.w_dvdp_cross1 
    mom_cross2 = data.u_dvdx_cross2 + data.v_dvdy_cross2 + data.w_dvdp_cross2
    mom_trans = data.uv_trans_dx + data.vv_trans_dy + data.vw_trans_dp
    mom_stat = data.u_dvdx_stat + data.v_dvdy_stat + data.w_dvdp_stat
    
    mom_crossu  = data.u_dvdx_cross1 + data.u_dvdx_cross2
    mom_crossv  = data.v_dvdy_cross1 + data.v_dvdy_cross2
    mom_crossw  = data.w_dvdp_cross1 + data.w_dvdp_cross2
    
        
    mom_sum = fu_local + fu_mean + dphidy + mom_mean + mom_trans + mom_stat + mom_cross + metric
    
    #levels = np.arange(-40,41.1,4.)
    levels = np.arange(-300.,300.1,30.)
    
    mn_dic = month_dic(1)
    tickspace = list(range(13,72,18))
    ticklabels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    # Nine subplots
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    
    plot_vars = [fu_mean, mom_mean, mom_sum, mom_trans, dphidy, metric]
    axes = [ax1,ax2,ax3,ax4,ax5,ax6]
    labels = ['a)','b)','c)','d)','e)','f)']
    for i in range(len(plot_vars)):
        f1=plot_vars[i].plot.contourf(ax=axes[i], x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
        if plot_precip:
            totp.plot.contour(ax=axes[i], x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
        axes[i].set_ylim(-60,60)
        axes[i].grid(True,linestyle=':')
        axes[i].set_yticks(np.arange(-60.,61.,30.))
        axes[i].set_xticks(tickspace)
    
    for i in [0,3]:
        axes[i].set_ylabel('Latitude')
        axes[i].text(-15, 60, labels[i])
    
    for i in [1,2,4,5]:
        axes[i].text(-5, 60, labels[i])
    
    for i in [3,4,5]:
        axes[i].set_xticklabels(ticklabels,rotation=25)
        
    # set titles    
    ax1.set_title('Zonal mean Coriolis', fontsize=17)
    ax2.set_title('Mean state advection', fontsize=17)
    ax3.set_title('Residual', fontsize=17)
    ax4.set_title('Transient eddy flux conv.', fontsize=17)
    ax5.set_title('Geopotential gradient', fontsize=17)
    ax6.set_title('Metric term', fontsize=17)
    
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.25, wspace=0.12)
    #Colorbar
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)

    
    if lonin == [-1.,361.]:
        figname = 'merid_mom_budg_' +run+ '.pdf'
    else:
        figname = 'merid_mom_budg_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()



mom_budg_hm('ap_2')
mom_budg_hm('ap_20')
mom_budg_hm('half_shallow')
mom_budg_hm('half_shallow', lonin=[90.,180.])
mom_budg_hm('half_shallow', lonin=[270.,360.])




