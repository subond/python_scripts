# Evaluate and plot the momentum budget terms over a given latitude range

from data_handling import time_means, month_dic, pentad_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from finite_difference import cfd
import gc

plt.rc('text', usetex=True)
font = {'size'   : 18}
plt.rc('font', **font)
mathtxt={'fontset':'custom', 'default':'regular'}
plt.rc('mathtext',**mathtxt)

#data_p = time_means('q_flux_test', [121,481], filename='plev_pentad', timeav='pentad')
data = xr.open_dataset( 'q_flux_test.nc', decode_times=False)
#data['totp'] = (('xofyear','lat','lon'), (data_p.convection_rain + data_p.condensation_rain)*86400.)
#data.totp.load()
#data_p.close()
land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/topo_10m/input/land.nc'
land = xr.open_dataset( land_file)
land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
        
a= 6376.0e3 #radius used in model
coslat = np.cos(data.lat * np.pi/180)    
omega = 7.2921150e-5
f = 2 * omega * np.sin(data.lat *np.pi/180)
mn_dic = month_dic(1)
pd_dic = pentad_dic(1)

print 'ready to start!'

def u_gradients(data):
        
    data['dudx'] = (('xofyear','pfull','lat','lon'),  cfd(data.ucomp.values,data.lon*np.pi/180.,3) )
    data['dudx'] = (('xofyear','pfull','lat','lon'), data.dudx/coslat/a )
    
    data['dudy'] = (('xofyear','pfull','lat','lon'),  cfd((data.ucomp*coslat).values,data.lat*np.pi/180.,2) )
    data['dudy'] = (('xofyear','pfull','lat','lon'), data.dudy/coslat/a )
    
    data['dudp'] = (('xofyear','pfull','lat','lon'),  cfd(data.ucomp.values,data.pfull*100.,1) )

def u_transients(data,lons):
    
    uu_ed = data.ucomp_sq - data.ucomp*data.ucomp
    uv_ed = data.ucomp_vcomp - data.ucomp*data.vcomp
    uw_ed = data.ucomp_omega - data.ucomp*data.omega
    
    data['duudx_ed'] = (('xofyear','pfull','lat','lon'),  cfd(uu_ed.values,data.lon*np.pi/180.,3) )
    data['duudx_ed'] = (('xofyear','pfull','lat','lon'), data.duudx_ed/coslat/a )
    
    data['duvdy_ed'] = (('xofyear','pfull','lat','lon'),  cfd((uv_ed*coslat*coslat).values,data.lat*np.pi/180.,2) )
    data['duvdy_ed'] = (('xofyear','pfull','lat','lon'), data.duvdy_ed/coslat/coslat/a )
    
    data['duwdp_ed'] = (('xofyear','pfull','lat','lon'),  cfd(uw_ed.values,data.pfull*100.,1) )
    
    data['mom_trans'] = (('xofyear','pfull','lat'), -1.*(data.duudx_ed[:,:,:,lons] + data.duvdy_ed[:,:,:,lons] + data.duwdp_ed[:,:,:,lons]).mean('lon') )
    
def zed(a,lons):
    a_zmean = a[:,:,:,lons].mean(('lon'))
    a_zed = a[:,:,:,lons] - a_zmean
    return a_zed


def mean_state_adv(data,lons):
    
    data['mom_mean'] = (('xofyear','pfull','lat'), -data.ucomp[:,:,:,lons].mean('lon') * data.dudx[:,:,:,lons].mean('lon') 
                       -  data.vcomp[:,:,:,lons].mean('lon') * data.dudy[:,:,:,lons].mean('lon') 
                       -  data.omega[:,:,:,lons].mean('lon') * data.dudp[:,:,:,lons].mean('lon') )

def stat_adv(data,lons):
    
    data['mom_stat'] = (('xofyear','pfull','lat'), -1.*(zed(data.ucomp,lons) * zed(data.dudx,lons) 
                       +  zed(data.vcomp,lons) * zed(data.dudy,lons) 
                       +  zed(data.omega,lons) * zed(data.dudp,lons)).mean('lon') )
                       
                       
def fv_and_dphidx(data,lons):
    data['dphidx'] = (('xofyear','pfull','lat','lon'),  cfd(data.height.values,data.lon*np.pi/180.,3) )
    data['dphidx'] = (('xofyear','pfull','lat'), data.dphidx[:,:,:,lons].mean('lon')/coslat/a*-9.8 )
    
    data['fv'] = (('xofyear','pfull','lat'), data.vcomp[:,:,:,lons].mean('lon')*f )
    
    
def plot_mom_var(data,var,lev,levels,lons):

    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]                
    plot_data = data.data_vars[var]
    plot_data = plot_data*10000.
    g=plot_data[:,lev,:].plot.contourf(x='xofyear', y='lat',levels=levels, add_label = False, add_colorbar=False, extend='both')
    cb1=plt.colorbar(g)
    cb1.set_label('$\displaystyle10^{-4}m/s^2$')
    cs=data.totp[:,:,lons].mean('lon').plot.contour(x='xofyear', y='lat',levels=np.arange(6.,31.,6.), add_label = False, add_colorbar=False,colors='k')
    plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
    plt.ylim(-45,45)
    plt.xlim(1,73)
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.ylabel('Latitude')
    plt.tight_layout()
    plt.title('')
    plt.savefig('/scratch/rg419/plots/qflux_run_analysis/topo/'+var+'_lattime.png')
    plt.close()




schneider = [i for i in range(len(data.lon)) if data.lon[i] >= 70. and data.lon[i] < 100.]
africa = [i for i in range(len(data.lon)) if data.lon[i] >= 330. or data.lon[i] < 60.]
asia = [i for i in range(len(data.lon)) if data.lon[i] >= 60. and data.lon[i] < 150.]
pacific = [i for i in range(len(data.lon)) if data.lon[i] >= 150. and data.lon[i] < 240.]
america = [i for i in range(len(data.lon)) if data.lon[i] >= 240. and data.lon[i] < 330.]

lons=asia

u_gradients(data)
print 'gradients done'
u_transients(data,lons)
print 'transients done'
mean_state_adv(data,lons)
print 'mean done'
stat_adv(data,lons)
print 'stat done'
fv_and_dphidx(data,lons)
print 'fv done, starting plotting'

tickspace = range(13,72,18)
labels = [mn_dic[(k+5)/6 ] for k in tickspace]

vars = ['fv','mom_trans','mom_stat','mom_mean','dphidx']
for var in vars:
    plot_mom_latlon(data,var,9,np.arange(-2.,2.1,0.2))
    #plot_mom_var(data,var,9,np.arange(-2.,2.1,0.2),lons)   
        
#mom_sum = (data.fv + data.dphidx + data.mom_trans+ data.mom_mean + data.mom_stat)*10000. 
#mom_sum[:,9,:].plot.contourf(x='xofyear', y='lat', levels = np.arange(-2.,2.1,0.2))
