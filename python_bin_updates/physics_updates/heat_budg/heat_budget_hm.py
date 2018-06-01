'''29/05/2018 Plot up terms contributing to heat budget'''

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling_updates import month_dic, gradients as gr
import sh
from pylab import rcParams
from climatology import precip_centroid
from hadley_cell import mass_streamfunction

def partition_advection(data, lons, lev=850, do_theta=False):
    
    if do_theta:
        convTtotheta=(1000./data.pfull)**(2./7.)
        theta = data.temp*convTtotheta
        data['ucomp_temp'] = data['ucomp_temp']*convTtotheta
        data['vcomp_temp'] = data['vcomp_temp']*convTtotheta
        data['omega_temp'] = data['omega_temp']*convTtotheta
        data['temp'] = data['temp']*convTtotheta
    
    #First do uT terms
    uT_trans_dx = -86400. * gr.ddx( (data.ucomp_temp - data.ucomp*data.temp).sel(pfull=lev) ) # <u'T'> = <uT> - <u><T>
    
    u = data.ucomp.sel(pfull=lev) # u
    T_dx = -86400. * gr.ddx( data.temp.sel(pfull=lev) )  # dTdx
    
    u_ed = u - u.mean('lon')
    T_dx_ed = T_dx - T_dx.mean('lon')
    
    u_dTdx_zav = u.mean('lon') * T_dx.mean('lon') # [u][dTdx], where brackets denote mean over all longitudes
    
    u_dTdx_cross1 = (u.mean('lon') * T_dx_ed).sel(lon=lons).mean('lon') # [u]dTdx*
    
    u_dTdx_cross2 = (u_ed * T_dx.mean('lon')).sel(lon=lons).mean('lon') # u*[dTdx]

    u_dTdx_stat = (u_ed * T_dx_ed).sel(lon=lons).mean('lon')          # u*dTdx* 
    
    data['uT_trans_dx'] = (('xofyear','lat'), uT_trans_dx.sel(lon=lons).mean('lon') )	
    data['u_dTdx_cross1'] = (('xofyear','lat'), u_dTdx_cross1 )	
    data['u_dTdx_cross2'] = (('xofyear','lat'), u_dTdx_cross2 )	
    data['u_dTdx_stat'] = (('xofyear','lat'), u_dTdx_stat )	
    data['u_dTdx_zav']  = (('xofyear','lat'), u_dTdx_zav )
    
    print('uT terms done')
    
    #Next do vT terms
    vT_trans_dy = -86400. * gr.ddy( (data.vcomp_temp - data.vcomp*data.temp).sel(pfull=lev))

    v = data.vcomp.sel(pfull=lev).load() # v
    T_dy = -86400. * gr.ddy( data.temp.sel(pfull=lev), vector=False)  # dTdy
    
    v_ed = v - v.mean('lon')
    T_dy_ed = T_dy - T_dy.mean('lon')
    
    v_dTdy_zav = v.mean('lon') * T_dy.mean('lon') # [v][dTdy]
    
    v_dTdy_cross1 = (v.mean('lon') * T_dy_ed).sel(lon=lons).mean('lon') # [v]dTdy*
    
    v_dTdy_cross2 = (v_ed * T_dy.mean('lon')).sel(lon=lons).mean('lon') # v*[dTdy]
        
    v_dTdy_stat = (v_ed * T_dy_ed).sel(lon=lons).mean('lon')           # v*dTdy* 
        
    data['vT_trans_dy'] = (('xofyear','lat'), vT_trans_dy.sel(lon=lons).mean('lon'))	
    data['v_dTdy_cross1'] = (('xofyear','lat'), v_dTdy_cross1 )	
    data['v_dTdy_cross2'] = (('xofyear','lat'), v_dTdy_cross2 )
    data['v_dTdy_stat'] = (('xofyear','lat'), v_dTdy_stat)	
    data['v_dTdy_zav']  = (('xofyear','lat'), v_dTdy_zav )
    
    print('vT terms done')
    
    #Finally do uw terms
    wT_trans_dp = -86400. * gr.ddp( (data.omega_temp - data.omega * data.temp).sel(lon=lons).mean('lon') )
    
    w = data.omega.sel(pfull=lev).load() # w
    T_dp = -86400. * (gr.ddp(data.temp)).sel(pfull=lev)  # dTdp
    
    w_ed = w - w.mean('lon')
    T_dp_ed = T_dp - T_dp.mean('lon')
    
    w_dTdp_zav = w.mean('lon') * T_dp.mean('lon') # [w][dTdp]
    
    w_dTdp_cross1 = (w.mean('lon') * T_dp_ed).sel(lon=lons).mean('lon') # [w]dTdp*
    
    w_dTdp_cross2 = (w_ed * T_dp.mean('lon')).sel(lon=lons).mean('lon') # w*[dTdp]
    
    w_dTdp_stat = (w_ed * T_dp_ed).sel(lon=lons).mean('lon')         # w*dTdp* 
    
    data['wT_trans_dp'] = (('xofyear','lat'), wT_trans_dp.sel(pfull=lev))	
    data['w_dTdp_cross1'] = (('xofyear','lat'), w_dTdp_cross1 )	
    data['w_dTdp_cross2'] = (('xofyear','lat'), w_dTdp_cross2 )
    data['w_dTdp_stat'] = (('xofyear','lat'), w_dTdp_stat)	
    data['w_dTdp_zav']  = (('xofyear','lat'), w_dTdp_zav )	
    
    print('wT terms done')
    

    
    
def heat_budg_hm(run, lev=850, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.], plot_precip=True, do_theta=False):
    
    rcParams['figure.figsize'] = 12, 8
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/heat_budg/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    #advective terms
    partition_advection(data, lons, lev=lev)
    
    try:
        precip_centroid(data)
    except:
        data['precipitation'] = data.convection_rain + data.condensation_rain
        precip_centroid(data)
        
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
     
    heat_mean = data.u_dTdx_zav #+ data.w_dTdp_zav #+data.v_dTdy_zav #
    heat_cross = data.u_dTdx_cross1 + data.v_dTdy_cross1 + data.w_dTdp_cross1 + data.u_dTdx_cross2 + data.v_dTdy_cross2 + data.w_dTdp_cross2
    heat_cross1 = data.u_dTdx_cross1 + data.v_dTdy_cross1 + data.w_dTdp_cross1 
    heat_cross2 = data.u_dTdx_cross2 + data.v_dTdy_cross2 + data.w_dTdp_cross2
    heat_trans = data.uT_trans_dx + data.vT_trans_dy + data.wT_trans_dp
    heat_stat = data.u_dTdx_stat + data.v_dTdy_stat + data.w_dTdp_stat
    
    heat_crossu  = data.u_dTdx_cross1 + data.u_dTdx_cross2
    heat_crossv  = data.v_dTdy_cross1 + data.v_dTdy_cross2
    heat_crossw  = data.w_dTdp_cross1 + data.w_dTdp_cross2
    
    tdt_rad = data.tdt_rad.sel(pfull=lev).sel(lon=lons).mean('lon')*86400.
    tdt_conv = data.dt_tg_convection.sel(pfull=lev).sel(lon=lons).mean('lon')*86400.
    tdt_cond = data.dt_tg_condensation.sel(pfull=lev).sel(lon=lons).mean('lon')*86400.
    tdt_diff = data.dt_tg_diffusion.sel(pfull=lev).sel(lon=lons).mean('lon')*86400.
    
    if do_theta:
        convTtotheta=(1000./data.pfull)**(2./7.)
        diabatic = (tdt_rad + tdt_conv + tdt_cond + tdt_diff)*convTtotheta.sel(pfull=lev)
        heat_sum = heat_mean + heat_trans + heat_stat + heat_cross + diabatic
        asc_cool = heat_sum * 0.
        Tdt = gr.ddt(data.temp.sel(pfull=lev)).sel(lon=lons).mean('lon') *86400.
    else:
        diabatic = tdt_rad + tdt_conv + tdt_cond + tdt_diff
        rho = data.pfull*100./data.temp/287.04
        asc_cool = (data.omega /(1004.64 * rho)).sel(pfull=lev).sel(lon=lons).mean('lon') *86400.
        heat_sum = heat_mean + heat_trans + heat_stat + heat_cross + diabatic + asc_cool
        Tdt = gr.ddt(data.temp.sel(pfull=lev)).sel(lon=lons).mean('lon') *86400.
    
    #levels = np.arange(-20,21.1,2.)
    #levels = np.arange(-10,10.1,1.)
    levels = np.arange(-1,1.1,0.1)
    
    mn_dic = month_dic(1)
    tickspace = list(range(13,72,18))
    ticklabels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    # Nine subplots
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    
    plot_vars = [heat_mean, heat_trans, heat_sum, asc_cool, diabatic, 
                 tdt_rad, tdt_diff, tdt_conv, tdt_cond]
                 
    axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
    labels = ['a)','b)','c)','d)','e)','f)','g)','h)','i)']
    for i in range(9):
        f1=plot_vars[i].plot.contourf(ax=axes[i], x='xofyear', y='lat', extend='both', add_labels=False, cmap='RdBu_r', levels = levels, add_colorbar=False)
        if plot_precip:
            psi.sel(pfull=500).plot.contour(ax=axes[i], x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
            psi.sel(pfull=500).plot.contour(ax=axes[i], x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
            psi.sel(pfull=500).plot.contour(ax=axes[i], x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
            data.p_cent.plot.line(ax=axes[i],color='k', linewidth=2)
            axes[i].set_xlabel(' ')
            axes[i].set_ylabel(' ')
            #totp.plot.contour(ax=axes[i], x='xofyear', y='lat', extend='both', levels=np.arange(-92.,109.,100.), add_colorbar=False, add_labels=False, alpha=0.25, colors='k', linewidths=2)
        axes[i].set_ylim(-60,60)
        axes[i].grid(True,linestyle=':')
        axes[i].set_yticks(np.arange(-60.,61.,30.))
        axes[i].set_xticks(tickspace)
    
    for i in [0,3,6]:
        axes[i].set_ylabel('Latitude')
        axes[i].text(-15, 60, labels[i])
    
    for i in [1,2,4,5,7,8]:
        axes[i].text(-5, 60, labels[i])
    
    for i in [6,7,8]:
        axes[i].set_xticklabels(ticklabels,rotation=25)
        
    # set titles    
    ax1.set_title('Mean state advection', fontsize=17)
    ax2.set_title('Transient eddy flux conv.', fontsize=17)
    ax3.set_title('Residual', fontsize=17)
    ax4.set_title('w cooling', fontsize=17)
    ax5.set_title('Total diabatic heating', fontsize=17)
    ax6.set_title('Radiative heating', fontsize=17)
    ax7.set_title('Diffusive heating', fontsize=17)
    ax8.set_title('Convective heating', fontsize=17)
    ax9.set_title('Condensational heating', fontsize=17)
    
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.25, wspace=0.12)
    #Colorbar
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)

    
    if lonin == [-1.,361.]:
        figname = 'zon_heat_budg_' +run+ '_u_850.pdf'
    else:
        figname = 'zon_heat_budg_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
   
if __name__ == "__main__":
    
    heat_budg_hm('ap_2',lev=850)

