"""
Evaluate and plot momentum budget at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means
from physics import mass_streamfunction
import sh
from finite_difference import cfd
from pylab import rcParams


def partitioning(field1, field2, fieldprod, lonmean=True):
    field1_zav = field1.mean('lon')
    field2_zav = field2.mean('lon')
    
    field12_trans = fieldprod - field1 * field2
    field12_stat = field1 * field2 - field1_zav * field2_zav
    
    if lonmean:
        field12_trans = field12_trans.mean('lon')
        field12_stat  = field12_stat.mean('lon')
    
    return field12_trans, field12_stat


def ddx(field):
    a= 6376.0e3 #radius used in model
    coslat = np.cos(field.lat * np.pi/180)
   
    field_dx = xr.DataArray( cfd( field.values, field.lon*np.pi/180, 2 ), [('xofyear', field.xofyear ), ('lat', field.lat), ('lon', field.lon)] )
    field_dx = -86400.*field_dx/coslat/a
    return field_dx
    
def ddy(field, prod=True):
    a= 6376.0e3 #radius used in model
    coslat = np.cos(field.lat * np.pi/180)
    cosfac = coslat**2. if prod==True else coslat
    
    field_cosfac = (field*cosfac)
    
    field_dy = xr.DataArray( cfd( field_cosfac.values, field.lat*np.pi/180, 1 ), [('xofyear', field.xofyear ), ('lat', field.lat)])
    field_dy = -86400.*field_dy/cosfac/a
    return field_dy
    
def ddp(field):
    field_dp = xr.DataArray( cfd( field.values, field.pfull*100., 1 ), [('xofyear', field.xofyear ), ('pfull', field.pfull), ('lat', field.lat)])
    field_dp = -86400.*field_dp
    return field_dp


def flux_conv(data, lons, lev=17):
    
    uu_trans, uu_stat = partitioning(data.ucomp[:,:,:,lons], data.ucomp[:,:,:,lons], data.ucomp_sq[:,:,:,lons], lonmean=False)
    uv_trans, uv_stat = partitioning(data.ucomp[:,:,:,lons], data.vcomp[:,:,:,lons], data.ucomp_vcomp[:,:,:,lons])
    uw_trans, uw_stat = partitioning(data.ucomp[:,:,:,lons], data.omega[:,:,:,lons], data.ucomp_omega[:,:,:,lons])
    
    print uu_stat.get_axis_num('lon')
    
    uu_trans_dx = ddx(uu_trans[:,lev,:,:])
    uv_trans_dy = ddy(uv_trans[:,lev,:])
    uw_trans_dp = ddp(uw_trans)

    uu_stat_dx = ddx(uu_stat[:,lev,:,:])
    uv_stat_dy = ddy(uv_stat[:,lev,:])
    uw_stat_dp = ddp(uw_stat)
    
    u_lev = data.ucomp[:,lev,:,:].load()
    v_lev = data.vcomp[:,lev,:,:].load()
    w_lev = data.omega[:,lev,:,:].load()
    
    u_dx = ddx(u_lev)
    u_dy = ddy(u_lev[:,:,lons].mean('lon'), prod=False)
    u_dp = ddp(data.ucomp[:,:,:,lons].mean('lon'))
    
    data['uu_trans_dx'] = (('xofyear','lat'), uu_trans_dx.mean('lon'))	
    data['uu_stat_dx'] = (('xofyear','lat'), uu_stat_dx.mean('lon'))	
    data['u_dudx'] = (('xofyear','lat'), (u_lev[:,:,lons]*u_dx[:,:,lons]).mean('lon') )	
    
    data['uv_trans_dy'] = (('xofyear','lat'), uv_trans_dy)	
    data['uv_stat_dy'] = (('xofyear','lat'), uv_stat_dy)	
    data['v_dudy'] = (('xofyear','lat'), v_lev[:,:,lons].mean('lon')*u_dy )	
    
    data['uw_trans_dp'] = (('xofyear','lat'), uw_trans_dp[:,lev,:])	
    data['uw_stat_dp'] = (('xofyear','lat'), uw_stat_dp[:,lev,:])	
    data['w_dudp'] = (('xofyear','lat'), w_lev[:,:,lons].mean('lon')*u_dp[:,lev,:] )	
    
    
def mom_budg_hm(run, months, lev=17, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 15, 10
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)

    if lonin[1]>lonin[0]:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    #advective terms
    flux_conv(data, lons, lev=lev)
    
    #Coriolis
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    fv = data.vcomp[:,lev,:,:] * f * 86400.
    fv = fv[:,:,lons].mean('lon')
    
    #Geopotential gradient
    dphidx = ddx(data.height[:,lev,:,:])
    dphidx = 9.8*dphidx[:,:,lons].mean('lon')
        
    mom_mean = data.u_dudx + data.v_dudy + data.w_dudp 
    mom_trans = data.uu_trans_dx + data.uv_trans_dy + data.uw_trans_dp
    mom_stat = data.uu_stat_dx + data.uv_stat_dy + data.uw_stat_dp
    
    mom_sum = fv + dphidx + mom_mean + mom_trans + mom_stat
    
    levels = np.arange(-10,10.1,2.)
    
    # Six subplots
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
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
    ax4.set_xlabel('Pentad')
    ax4.set_ylim(-60,60)
    
    #Fifth plot
    mom_stat.plot.contourf(ax=ax5, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax5.grid(True,linestyle=':')
    ax5.set_xlabel('Pentad')
    ax5.set_ylim(-60,60)
    
    #Sixth plot
    mom_sum.plot.contourf(ax=ax6, x='xofyear', y='lat', extend='both', levels = levels, add_colorbar=False, add_labels=False)
    ax6.grid(True,linestyle=':')
    ax6.set_xlabel('Pentad')
    ax6.set_ylim(-60,60)
    
    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
    #Colorbar
    cb1=f.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
    cb1.set_label('$ms^{-1}day^{-1}$')
    
    plt.savefig(plot_dir+'zon_mom_budg_60_150.pdf', format='pdf')
    plt.close()
        


#mom_budg_hm('ap_2', [121,481])
mom_budg_hm('full_qflux', [121,133])#, lonin=[60.,150.])
#mom_budg_hm('flat_qflux', [121,481], lonin=[60.,150.])


