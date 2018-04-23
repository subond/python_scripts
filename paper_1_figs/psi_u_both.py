"""
Function to plot meridional overturning, zonal wind speed, and eddy flux convergences

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means
from physics import mass_streamfunction
import sh
from finite_difference import cfd_old as cfd
from pylab import rcParams

    
rcParams['figure.figsize'] = 12, 7
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def uv_partitioning(data, start_index, lons):
    
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    
    u_ztav = data.ucomp[start_index:start_index+4,:,:,:].sel(lon=lons).mean(('xofyear','lon')) 
    v_ztav = data.vcomp[start_index:start_index+4,:,:,:].sel(lon=lons).mean(('xofyear','lon'))
    
    uv_eddy = data.ucomp_vcomp[start_index:start_index+4,:,:,:].sel(lon=lons).mean(('xofyear','lon')) - u_ztav * v_ztav 
                     
    uv_conv = xr.DataArray( cfd( (uv_eddy*coslat*coslat).values, data.lat*np.pi/180, 1 ), [('pfull', data.pfull ), ('lat', data.lat)])
    uv_conv = -86400.*uv_conv/coslat/coslat/a

    return uv_conv, u_ztav
    
    
def u_psi(run, bef_aft, lonin=[-1.,361.]):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    print data.vcomp.sel(pfull=900.).sel(xofyear=40.)[50,33].values
    vtest = np.nan_to_num(data.vcomp)
    vtest = xr.DataArray(vtest, [('xofyear', data.xofyear), ('pfull', data.pfull), ('lat', data.lat), ('lon', data.lon)])
    
    #data['vcomp'] = (('xofyear','pfull','lat','lon'), vtest)
    print data.vcomp.sel(pfull=900.).sel(xofyear=40.)[50,33].values

    
    psi = mass_streamfunction(data, a=6376.0e3, lons=lons, dp_in=50.)
    psi /= 1.e9
    
    psi_before = psi[:,bef_aft[0]:bef_aft[0]+4,::-1].mean('xofyear').transpose()
    psi_after = psi[:,bef_aft[1]:bef_aft[1]+4,::-1].mean('xofyear').transpose()
    
    uv_conv_before, u_ztav_before = uv_partitioning(data, bef_aft[0], lons)
    uv_conv_after, u_ztav_after  = uv_partitioning(data, bef_aft[1], lons)
    
    ds = xr.Dataset({'psi_before': (['pfull', 'lat'], psi_before),
                     'psi_after':  (['pfull', 'lat'], psi_after),
                     'uv_conv_before':  (['pfull', 'lat'], uv_conv_before),
                     'uv_conv_after':  (['pfull', 'lat'], uv_conv_after),
                     'u_ztav_before':  (['pfull', 'lat'], u_ztav_before),
                     'u_ztav_after':   (['pfull', 'lat'], u_ztav_after)},
                     coords={'pfull': ('pfull', data.pfull),
                               'lat': ('lat', data.lat)})
    
    return ds


data_ap2 = u_psi('ap_2', [30,39])
data_full = u_psi('full_qflux', [18,39], [60.,150.])


lwid=2


# Four subplots
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
plt.set_cmap('RdBu_r')
#First plot
#f1 = ax1.contourf(data_ap2.lat, data_ap2.pfull, data_ap2.uv_conv_before, extend = 'both', levels = np.arange(-10.,10.1,1.))
f1 = data_ap2.u_ztav_before.plot.contourf(x='lat', y='pfull', ax=ax1, levels=np.arange(-60.,61.,10.), extend='both', add_labels=False, add_colorbar=False, yincrease=False)
#ax1.contourf(data_ap2.lat, data_ap2.pfull, data_ap2.u_ztav_before, levels=np.arange(-60.,61.,10.))
ax1.contour(data_ap2.lat, data_ap2.pfull, data_ap2.psi_before, colors='k', levels=np.arange(-300.,301.,60.), linewidths=lwid)
#ax1.invert_yaxis()
ax1.set_title('$ap2$, pre-onset', fontsize=17)
ax1.set_ylabel('Pressure, hPa')
ax1.set_xlim(-60,60)
ax1.grid(True,linestyle=':')
ax1.text(-80, 0, 'a)')

#Second plot
#f2 = ax2.contourf(data_ap2.lat, data_ap2.pfull, data_ap2.uv_conv_after, extend = 'both', levels = np.arange(-10.,10.1,1.))
f2 = data_ap2.u_ztav_after.plot.contourf(x='lat', y='pfull', ax=ax2, levels=np.arange(-60.,61.,10.), extend='both', add_labels=False, add_colorbar=False, yincrease=False)
#f2 = ax2.contourf(data_ap2.lat, data_ap2.pfull, data_ap2.u_ztav_after, levels=np.arange(-60.,61.,10.))
ax2.contour(data_ap2.lat, data_ap2.pfull, data_ap2.psi_after, colors='k', levels=np.arange(-300.,301.,60.), linewidths=lwid)
ax2.set_xlim(-60,60)
ax2.set_title('$ap2$, post-onset', fontsize=17)
ax2.grid(True,linestyle=':')
ax2.text(-70, 0, 'b)')

#Third plot
#f3 = ax3.contourf(data_full.lat, data_full.pfull, data_full.uv_conv_before, extend = 'both', levels = np.arange(-10.,10.1,1.))
f3 = data_full.u_ztav_before.plot.contourf(x='lat', y='pfull', ax=ax3, levels=np.arange(-60.,61.,10.), extend='both', add_labels=False, add_colorbar=False, yincrease=False)
#f3 = ax3.contourf(data_full.lat, data_full.pfull, data_full.u_ztav_before, levels=np.arange(-60.,61.,10.))
ax3.contour(data_full.lat, data_full.pfull, data_full.psi_before, colors='k', levels=np.arange(-300.,301.,60.), linewidths=lwid)
ax3.grid(True,linestyle=':')
#ax3.invert_yaxis()
ax3.set_title('$full$, pre-onset', fontsize=17)
ax3.set_ylabel('Pressure, hPa')
ax3.set_xlabel('Latitude')
ax3.text(-80, 0, 'c)') 

#Fourth plot
#f2 = ax4.contourf(data_full.lat, data_full.pfull, data_full.uv_conv_after, extend = 'both', levels = np.arange(-10.,10.1,1.))
f4 = data_full.u_ztav_after.plot.contourf(x='lat', y='pfull', ax=ax4, levels=np.arange(-60.,61.,10.), extend='both', add_labels=False, add_colorbar=False, yincrease=False)
#f4 = ax4.contourf(data_full.lat, data_full.pfull, data_full.u_ztav_after, levels=np.arange(-60.,61.,10.))
ax4.contour(data_full.lat, data_full.pfull, data_full.psi_after, colors='k', levels=np.arange(-300.,301.,60.), linewidths=lwid)
ax4.set_xlim(-60,60)
ax4.set_title('$full$, post-onset', fontsize=17)
ax4.grid(True,linestyle=':')
ax4.set_xlabel('Latitude')
ax4.text(-70, 0, 'd)')

    
#plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=f.colorbar(f2, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.12, aspect=30, shrink=0.5)

plt.savefig(plot_dir+'psi_u_both.pdf', format='pdf')
plt.close()
