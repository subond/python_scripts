"""
Function to plot meridional overturning, zonal wind speed, and eddy flux convergences

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means
from physics import mass_streamfunction
import sh
from finite_difference import cfd
from pylab import rcParams

    
rcParams['figure.figsize'] = 12, 9
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
        
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/full_qflux.nc')
a= 6376.0e3 #radius used in model
coslat = np.cos(data.lat * np.pi/180)

lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= 60. and data.lon[i] < 150.]
  
psi = mass_streamfunction(data, a=6376.0e3, lons=lons, dp_in=50.)
psi /= 1.e9
    
psi_before = psi[::-1,18:18+4,::-1].mean('xofyear').transpose()
psi_after = psi[::-1,39:39+4,::-1].mean('xofyear').transpose()



def uv_partitioning(data, start_index):
    u_tav = data.ucomp[start_index:start_index+4,:,:,:].mean('xofyear')
    u_ztav = u_tav.sel(lon=lons).mean('lon')
    
    v_tav = data.vcomp[start_index:start_index+4,:,:,:].mean('xofyear')
    v_ztav = v_tav.sel(lon=lons).mean('lon')
    
    uv_trans = (data.ucomp_vcomp[start_index:start_index+4,:,:,:].mean('xofyear') 
                 - u_tav * v_tav ).sel(lon=lons).mean('lon')
                 
    uv_stat = (u_tav * v_tav).sel(lon=lons).mean('lon') - u_ztav * v_ztav
    
    uv_conv_trans = xr.DataArray( cfd( (uv_trans*coslat*coslat).values, data.lat*np.pi/180, 1 ), [('pfull', data.pfull ), ('lat', data.lat)])
    uv_conv_trans = -86400.*uv_conv_trans/coslat/coslat/a
    
    uv_conv_stat = xr.DataArray( cfd( (uv_stat*coslat*coslat).values, data.lat*np.pi/180, 1 ), [('pfull', data.pfull ), ('lat', data.lat)])
    uv_conv_stat = -86400.*uv_conv_stat/coslat/coslat/a
    
    return uv_conv_trans, uv_conv_stat, u_ztav
        
uv_conv_trans_before, uv_conv_stat_before, u_ztav_before = uv_partitioning(data, 18)
uv_conv_trans_after,  uv_conv_stat_after, u_ztav_after  = uv_partitioning(data, 39)

lwid=2

# Four subplots
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
plt.set_cmap('RdBu_r')
#First plot
f1 = ax1.contourf(data.lat, data.pfull, uv_conv_trans_before, extend = 'both', levels = np.arange(-10.,10.1,1.))
ax1.contour(data.lat, data.pfull, u_ztav_before, colors='0.7', levels=np.arange(-60.,61.,10.), linewidths=lwid)
ax1.contour(data.lat, data.pfull, psi_before, colors='k', levels=np.arange(-300.,301.,60.), linewidths=lwid)
ax1.invert_yaxis()
ax1.set_ylabel('Pressure, hPa')
ax1.set_xlim(-60,60)
ax1.grid(True,linestyle=':')
ax1.text(-80, 0, 'a)')

#Second plot
f1 = ax2.contourf(data.lat, data.pfull, uv_conv_stat_before, extend = 'both', levels = np.arange(-10.,10.1,1.))
ax2.contour(data.lat, data.pfull, u_ztav_before, colors='0.7', levels=np.arange(-60.,61.,10.), linewidths=lwid)
ax2.contour(data.lat, data.pfull, psi_before, colors='k', levels=np.arange(-300.,301.,60.), linewidths=lwid)
#ax2.invert_yaxis()
ax2.set_xlim(-60,60)
ax2.grid(True,linestyle=':')
ax2.text(-80, 0, 'b)')

#Third plot
f2 = ax3.contourf(data.lat, data.pfull, uv_conv_trans_after, extend = 'both', levels = np.arange(-10.,10.1,2.))
ax3.contour(data.lat, data.pfull, u_ztav_after, colors='0.7', levels=np.arange(-60.,61.,10.), linewidths=lwid)
ax3.contour(data.lat, data.pfull, psi_after, colors='k', levels=np.arange(-300.,301.,60.), linewidths=lwid)
ax3.grid(True,linestyle=':')
ax3.invert_yaxis()
ax3.set_ylabel('Pressure, hPa')
ax3.set_xlabel('Latitude')
ax3.text(-80, 0, 'c)') 

#Fourth plot
f1 = ax4.contourf(data.lat, data.pfull, uv_conv_stat_after, extend = 'both', levels = np.arange(-10.,10.1,1.))
ax4.contour(data.lat, data.pfull, u_ztav_after, colors='0.7', levels=np.arange(-60.,61.,10.), linewidths=lwid)
ax4.contour(data.lat, data.pfull, psi_after, colors='k', levels=np.arange(-300.,301.,60.), linewidths=lwid)
#ax4.invert_yaxis()
ax4.set_xlim(-60,60)
ax4.grid(True,linestyle=':')
ax4.set_xlabel('Latitude')
ax4.text(-80, 0, 'd)')

    
#plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=f.colorbar(f2, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30, shrink=0.5)
cb1.set_label('Eddy momentum flux convergence, ms$^{-1}$day$^{-1}$')

plt.savefig(plot_dir+'psi_u_full.pdf', format='pdf')
plt.close()
    


