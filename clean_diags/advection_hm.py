"""
Evaluate and plot advective terms in momentum budget in terms of contributions from uu, uv, uw products

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
    
    
    
def advection_hm(run, months, lev=150, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 15, 10
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    #advective terms
    partition_advection(data, lons, lev=150)
        
    data['mom_mean']  = (('xofyear','lat'), data.u_dudx_zav + data.v_dudy_zav + data.w_dudp_zav )
    data['mom_trans'] = (('xofyear','lat'), data.uu_trans_dx + data.uv_trans_dy + data.uw_trans_dp )
    data['mom_stat']  = (('xofyear','lat'), data.u_dudx_stat + data.v_dudy_stat + data.w_dudp_stat )
    
    levels = np.arange(-10,10.1,2.)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    fig_dicts = [
        {'p1':'mom_mean', 'p2':'u_dudx_zav', 'p3':'v_dudy_zav', 'p4':'w_dudp_zav', 'namebase':'meanstate_adv'},
        {'p1':'mom_stat', 'p2':'u_dudx_stat', 'p3':'v_dudy_stat', 'p4':'w_dudp_stat', 'namebase':'stat_adv'},
        {'p1':'mom_trans', 'p2':'uu_trans_dx', 'p3':'uv_trans_dy', 'p4':'uw_trans_dp', 'namebase':'trans_adv'}
        ]
        
    for fig_dict in fig_dicts:
        # Six subplots
        fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
        plt.set_cmap('RdBu_r')
        #First plot
        f1=data[fig_dict['p1']].plot.contourf(ax=ax1, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
        ax1.set_ylabel('Latitude')
        ax1.set_ylim(-60,60)
        ax1.grid(True,linestyle=':')
        
        #Fourth plot
        data[fig_dict['p2']].plot.contourf(ax=ax4, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
        ax4.grid(True,linestyle=':')
        ax4.set_ylabel('Latitude')
        ax4.set_xticks(tickspace)
        ax4.set_xticklabels(labels,rotation=25)
        ax4.set_ylim(-60,60)
        
        #Fifth plot
        data[fig_dict['p3']].plot.contourf(ax=ax5, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
        ax5.grid(True,linestyle=':')
        ax5.set_xticks(tickspace)
        ax5.set_xticklabels(labels,rotation=25)
        ax5.set_ylim(-60,60)
        
        #Sixth plot
        data[fig_dict['p4']].plot.contourf(ax=ax6, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
        ax6.grid(True,linestyle=':')
        ax6.set_xticks(tickspace)
        ax6.set_xticklabels(labels,rotation=25)
        ax6.set_ylim(-60,60)
        
        plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
        #Colorbar
        cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
        cb1.set_label('$ms^{-1}day^{-1}$')
        
        if lonin == [-1.,361.]:
            figname = fig_dict['namebase'] + '.pdf'
        else:
            figname = fig_dict['namebase'] + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
        
        plt.savefig(plot_dir + figname, format='pdf')
        plt.close()
    
    
advection_hm('am_qflux', [121,481])
advection_hm('am_qflux', [121,481], lonin=[60.,150.])
advection_hm('ap_20_qflux', [121,481])
advection_hm('ap_20_qflux', [121,481], lonin=[60.,150.])

#advection_hm('ap_2', [121,481])
#advection_hm('full_qflux', [121,481])
#advection_hm('full_qflux', [121,481], lonin=[60.,150.])
#advection_hm('flat_qflux', [121,481])
#advection_hm('flat_qflux', [121,481], lonin=[60.,150.])


