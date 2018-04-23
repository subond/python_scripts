# Make contour plots of the terms in the momentum budget at 200hPa plus the wind vectors for the topography run over the onset period. 

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

def u_transients(data):
    
    uu_ed = data.ucomp_sq - data.ucomp*data.ucomp
    uv_ed = data.ucomp_vcomp - data.ucomp*data.vcomp
    uw_ed = data.ucomp_omega - data.ucomp*data.omega
    
    data['duudx_ed'] = (('xofyear','pfull','lat','lon'),  cfd(uu_ed.values,data.lon*np.pi/180.,3) )
    data['duudx_ed'] = (('xofyear','pfull','lat','lon'), data.duudx_ed/coslat/a )
    
    data['duvdy_ed'] = (('xofyear','pfull','lat','lon'),  cfd((uv_ed*coslat*coslat).values,data.lat*np.pi/180.,2) )
    data['duvdy_ed'] = (('xofyear','pfull','lat','lon'), data.duvdy_ed/coslat/coslat/a )
    
    data['duwdp_ed'] = (('xofyear','pfull','lat','lon'),  cfd(uw_ed.values,data.pfull*100.,1) )
    
    data['mom_trans'] = (('xofyear','pfull','lat','lon'), -1.*(data.duudx_ed + data.duvdy_ed + data.duwdp_ed) )
    


def mean_state_adv(data):
    
    data['mom_mean'] = (('xofyear','pfull','lat','lon'), -data.ucomp * data.dudx 
                       -  data.vcomp * data.dudy 
                       -  data.omega * data.dudp )

                       
                       
def fv_and_dphidx(data):
    data['dphidx'] = (('xofyear','pfull','lat','lon'),  cfd(data.height.values,data.lon*np.pi/180.,3) )
    data['dphidx'] = (('xofyear','pfull','lat','lon'), data.dphidx/coslat/a*-9.8 )
    
    data['fv'] = (('xofyear','pfull','lat','lon'), data.vcomp*f )
    



def plot_mom_latlon(data,var,lev, levels):
        
    tickspace = range(65,150,20)
    
    plot_data = data.data_vars[var]
    plot_data = plot_data*10000.
    for mnth in range(0,12):
        g=plot_data[mnth*6 : mnth*6+6,lev,:,:].plot.contourf(x='lon', y='lat',levels=levels, col='xofyear', col_wrap=3, extend='both')
        plt.ylim(-45,45)
        plt.xlim(60,150)
        plt.xticks(tickspace)
        
        for i, ax in enumerate(g.axes.flat):
            land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
            ax.quiver(data.lon[::3], data.lat[::1], data.ucomp[mnth*6 +i,2,::1,::3], data.vcomp[mnth*6+i,2,::1,::3], scale=500.,headwidth=5)
            ax.set_title(pd_dic[mnth*6+i+1])
        plt.savefig('/scratch/rg419/plots/qflux_run_analysis/qflux/latlon/' + var + str(mnth + 1) + '.png')
        plt.close()
               

u_gradients(data)
print 'gradients done'
u_transients(data)
print 'transients done'
mean_state_adv(data)
print 'mean done'
fv_and_dphidx(data)
print 'fv done, starting plotting'

tickspace = range(13,72,18)
labels = [mn_dic[(k+5)/6 ] for k in tickspace]

vars = ['fv','mom_trans','mom_mean','dphidx']
for var in vars:
    plot_mom_latlon(data,var,9,np.arange(-4.,4.1,0.4))
