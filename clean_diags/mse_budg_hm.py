"""
Evaluate and plot moist static energy budget at 850 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams


def partition_mse_advection(data, lons, lev=850):
        
    cp = 287.04/2*7
    L = 2.500e6
    data['mse_u'] = cp*data.ucomp_temp + L*data.sphum_u
    data['mse_v'] = cp*data.vcomp_temp + L*data.sphum_v
    data['mse_w'] = cp*data.omega_temp + L*data.sphum_w
    data['mse'] = cp*data.temp + L*data.sphum
    
    #First do uu terms
    umse_trans_dx = -1. * gr.ddx( (data.mse_u - data.ucomp*data.mse).sel(pfull=lev) ) # <u'q'> = <uq> - <u><q>
    
    u = data.ucomp.sel(pfull=lev).load()
    mse = data.mse.sel(pfull=lev) # mse
    mse_dx = -1. * gr.ddx( mse )  # dmsedx
    
    u_dmsedx_zav = u.sel(lon=lons).mean('lon') * mse_dx.sel(lon=lons).mean('lon') # [u][dmsedx]

    u_dmsedx_stat = (u * mse_dx).sel(lon=lons).mean('lon') - u_dmsedx_zav           # u*dudx* = [udmsedx] - [u][dmsedx]
    
    data['umse_trans_dx'] = (('xofyear','lat'), umse_trans_dx.sel(lon=lons).mean('lon') )	
    data['u_dmsedx_stat'] = (('xofyear','lat'), u_dmsedx_stat )	
    data['u_dmsedx_zav']  = (('xofyear','lat'), u_dmsedx_zav )
    
    print 'umse terms done'
    
    #Next do uv terms
    vmse_trans_dy = -1. * gr.ddy( (data.mse_v - data.mse * data.vcomp).sel(pfull=lev) , uv=True)

    v = data.vcomp.sel(pfull=lev).load() # v
    mse_dy = -1. * gr.ddy( mse, vector=False)  # dmsedy
    
    v_dmsedy_zav = v.sel(lon=lons).mean('lon') * mse_dy.sel(lon=lons).mean('lon') # [v][dudy]
        
    v_dmsedy_stat = (v * mse_dy).sel(lon=lons).mean('lon') - v_dmsedy_zav           # v*dudy* = [vdudy] - [v][dudy]
        
    data['vmse_trans_dy'] = (('xofyear','lat'), vmse_trans_dy.sel(lon=lons).mean('lon'))	
    data['v_dmsedy_stat'] = (('xofyear','lat'), v_dmsedy_stat)	
    data['v_dmsedy_zav']  = (('xofyear','lat'), v_dmsedy_zav )
    
    print 'vmse terms done'
    
    #Finally do uw terms
    wmse_trans_dp = -1. * gr.ddp( (data.mse_w - data.mse * data.omega).sel(lon=lons).mean('lon') )
    
    w = data.omega.sel(pfull=lev).load() # w
    mse_dp = -1. * (gr.ddp(data.mse)).sel(pfull=lev)  # dudp
    
    w_dmsedp_zav = w.sel(lon=lons).mean('lon') * mse_dp.sel(lon=lons).mean('lon')
    w_dmsedp_stat = (w * mse_dp).sel(lon=lons).mean('lon') - w_dmsedp_zav
    
    data['wmse_trans_dp'] = (('xofyear','lat'), wmse_trans_dp.sel(pfull=lev))	
    data['w_dmsedp_stat'] = (('xofyear','lat'), w_dmsedp_stat)	
    data['w_dmsedp_zav']  = (('xofyear','lat'), w_dmsedp_zav )	
    
    print 'wmse terms done'
    
    
    
def mse_budg_hm(run, lev=850, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 15, 5
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
    
    #advective terms
    partition_mse_advection(data, lons, lev=150)
    
    cp = 287.04/2*7
    L = 2.500e6
    data['mse'] = (cp*data.temp + L*data.sphum).sel(lon=lons).mean('lon')
    mse_max_loc = [data.lat[i] for i in data.mse.sel(pfull=lev).argmax('lat').values]
    
    mse_mean = data.u_dmsedx_zav + data.v_dmsedy_zav + data.w_dmsedp_zav
    mse_trans = data.umse_trans_dx + data.vmse_trans_dy + data.wmse_trans_dp
    mse_stat = data.u_dmsedx_stat + data.v_dmsedy_stat + data.w_dmsedp_stat
    
    
    #levels = np.arange(-8000,8000.1,1000.)
    levels = np.arange(-0.16,0.161,0.02)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    # Three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
    plt.set_cmap('RdBu_r')
    #First plot
    f1=mse_mean.plot.contourf(ax=ax1, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax1.plot(data.xofyear,mse_max_loc, 'k')
    ax1.grid(True,linestyle=':')
    ax1.set_ylim(-60,60)
    ax1.set_ylabel('Latitude')
    ax1.set_xticks(tickspace)
    ax1.set_xticklabels(labels,rotation=25)    
    
    #Second plot
    mse_trans.plot.contourf(ax=ax2, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax2.plot(data.xofyear,mse_max_loc, 'k')
    ax2.grid(True,linestyle=':')
    ax2.set_xticks(tickspace)
    ax2.set_xticklabels(labels,rotation=25)
    ax2.set_ylim(-60,60)
    
    #Third plot
    mse_stat.plot.contourf(ax=ax3, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax3.plot(data.xofyear,mse_max_loc, 'k')
    ax3.grid(True,linestyle=':')
    ax3.set_xticks(tickspace)
    ax3.set_xticklabels(labels,rotation=25)
    ax3.set_ylim(-60,60)
    

    
    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
    #Colorbar
    cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30)
    cb1.set_label('$Jkg^{-1}day^{-1}$')
    
    if lonin == [-1.,361.]:
        figname = 'mse_budg.pdf'
    else:
        figname = 'mse_budg_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
        
#mse_budg_hm('ap_20_qflux')
#mse_budg_hm('ap_20_qflux', lonin=[60.,150.])
#mse_budg_hm('am_qflux')
#mse_budg_hm('am_qflux', lonin=[60.,150.])
mse_budg_hm('ap_2')
mse_budg_hm('full_qflux')
mse_budg_hm('full_qflux', lonin=[60.,150.])
#mom_budg_hm('flat_qflux', [121,481])
#mom_budg_hm('flat_qflux', [121,481], lonin=[60.,150.])


