"""
Evaluate and plot meridional momentum budget at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams

# Budget: Dv/Dt + fu + dphi/dy = F
# dv/dt = - udv/dx - vdv/dy - wdv/dp - fu - dphi/dy
# Zonal mean:
# dv/dt = - vdvdy - wdv/dp - fu - dphi/dy

def partition_advection(data, lons, lev=150):
        
    #Do vv terms
    vv_trans_dy = -86400. * gr.ddy( (data.vcomp_sq - data.vcomp * data.vcomp).sel(pfull=lev) , uv=True)

    v = data.vcomp.sel(pfull=lev).load() # v
    v_dy = -86400. * gr.ddy( v )  # dudy
    
    v_dvdy_zav = v.sel(lon=lons).mean('lon') * v_dy.sel(lon=lons).mean('lon') # [v][dudy]
        
    v_dvdy_stat = (v * v_dy).sel(lon=lons).mean('lon') - v_dvdy_zav           # v*dudy* = [vdudy] - [v][dudy]
        
    data['vv_trans_dy'] = (('xofyear','lat'), vv_trans_dy.sel(lon=lons).mean('lon'))	
    data['v_dvdy_stat'] = (('xofyear','lat'), v_dvdy_stat)	
    data['v_dvdy_zav']  = (('xofyear','lat'), v_dvdy_zav )
    
    print 'vv terms done'
    
    #Do vw terms
    vw_trans_dp = -86400. * gr.ddp( (data.vcomp_omega - data.vcomp * data.omega).sel(lon=lons).mean('lon') )
    
    w = data.omega.sel(pfull=lev).load() # w
    v_dp = -86400. * (gr.ddp(data.vcomp)).sel(pfull=lev)  # dudp
    
    w_dvdp_zav = w.sel(lon=lons).mean('lon') * v_dp.sel(lon=lons).mean('lon')
    w_dvdp_stat = (w * v_dp).sel(lon=lons).mean('lon') - w_dvdp_zav
    
    data['vw_trans_dp'] = (('xofyear','lat'), vw_trans_dp.sel(pfull=lev))	
    data['w_dvdp_stat'] = (('xofyear','lat'), w_dvdp_stat)	
    data['w_dvdp_zav']  = (('xofyear','lat'), w_dvdp_zav )	
    
    print 'vw terms done'
    
    
    
def mom_budg_hm(run, lev=150, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    a = 6376.0e3 
       
    rcParams['figure.figsize'] = 12, 5.5
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/itcz_theory_testing/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    #advective terms
    partition_advection(data, lons, lev=150)
    
    #Coriolis
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180.)
    fu = -1.*data.ucomp.sel(pfull=lev) * f * 86400.
    fu = fu.sel(lon=lons).mean('lon')
    
    metric_term = (-86400. * data.ucomp.sel(pfull=lev)*data.ucomp.sel(pfull=lev)*np.tan(data.lat*np.pi/180.)/a).mean('lon')
    
    #totp = ((data.convection_rain + data.condensation_rain)*86400.).sel(lon=lons).mean('lon')
    #abs_vort = (data.vor + f).sel(lon=lons).mean('lon')*86400.
    
    #Geopotential gradient
    dphidy = gr.ddy(data.height.sel(pfull=lev), vector=False)
    dphidy = -86400. * 9.8 * dphidy.sel(lon=lons).mean('lon')
        
    mom_mean = data.v_dvdy_zav + data.w_dvdp_zav
    mom_trans = data.vv_trans_dy + data.vw_trans_dp
    mom_stat =  data.v_dvdy_stat + data.w_dvdp_stat
    
    mom_sum = mom_mean + mom_trans + mom_stat + fu + dphidy + metric_term
    
    resid = -(fu + dphidy + metric_term)
    
    levels = np.arange(-10.,10.1,1.)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    resid.plot.contourf(x='xofyear', y='lat', extend='both', levels = levels, add_labels=False)
    plt.ylabel('Latitude')
    plt.ylim(-60,60)
    plt.xticks(tickspace,labels,rotation=25)
    plt.yticks(np.arange(-60.,61.,30.))
    plt.grid(True,linestyle=':')
    plt.tight_layout()  
    plt.savefig(plot_dir+'residual.pdf')
    plt.close()
    
    mom_mean.plot.contourf(x='xofyear', y='lat', extend='both', levels = levels, add_labels=False)
    plt.ylabel('Latitude')
    plt.ylim(-60,60)
    plt.xticks(tickspace,labels,rotation=25)
    plt.yticks(np.arange(-60.,61.,30.))
    plt.grid(True,linestyle=':')
    plt.tight_layout()  
    plt.savefig(plot_dir+'mom_mean.pdf')
    plt.close()


mom_budg_hm('ap_2')
