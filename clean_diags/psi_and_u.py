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


    
def upsi(run, months, before_after, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 15, 10
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    #data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
        
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    
    psi = mass_streamfunction(data, a=6376.0e3, lons=lons, dp_in=50.)
    psi /= 1.e9
    
    psi_before = psi[::-1,before_after[0]:before_after[0]+4,::-1].mean('xofyear').transpose()
    psi_after = psi[::-1,before_after[1]:before_after[1]+4,::-1].mean('xofyear').transpose()
    
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
        
    uv_conv_trans_before, uv_conv_stat_before, u_ztav_before = uv_partitioning(data, before_after[0])
    uv_conv_trans_after,  uv_conv_stat_after, u_ztav_after  = uv_partitioning(data, before_after[1])
    
    # Four subplots
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f1 = ax1.contourf(data.lat, data.pfull, uv_conv_trans_before, extend = 'both', levels = np.arange(-5.,5.1,1.))
    ax1.contour(data.lat, data.pfull, u_ztav_before, colors='0.75', levels=np.arange(-60.,61.,10.), linewidths=2)
    ax1.contour(data.lat, data.pfull, psi_before, colors='k', levels=np.arange(-300.,301.,50.), linewidths=2)
    ax1.invert_yaxis()
    ax1.set_ylabel('Pressure, hPa')
    ax1.set_xlim(-60,60)
    ax1.grid(True,linestyle=':')
    
    #Second plot
    f2 = ax2.contourf(data.lat, data.pfull, uv_conv_trans_after, extend = 'both', levels = np.arange(-5.,5.1,1.))
    ax2.contour(data.lat, data.pfull, u_ztav_after, colors='0.75', levels=np.arange(-60.,61.,10.))
    ax2.contour(data.lat, data.pfull, psi_after, colors='k', levels=np.arange(-300.,301.,50.), linewidths=2)
    ax2.grid(True,linestyle=':')
    ax2.set_xlim(-60,60)
    
    #Third plot
    f2 = ax3.contourf(data.lat, data.pfull, uv_conv_stat_before, extend = 'both', levels = np.arange(-5.,5.1,1.))
    ax3.contour(data.lat, data.pfull, u_ztav_before, colors='0.75', levels=np.arange(-60.,61.,10.), linewidths=2)
    ax3.contour(data.lat, data.pfull, psi_before, colors='k', levels=np.arange(-300.,301.,50.), linewidths=2)
    ax3.grid(True,linestyle=':')
    ax3.invert_yaxis()
    ax3.set_ylabel('Pressure, hPa')
    ax3.set_xlabel('Latitude')
    
    #Fourth plot
    f2 = ax4.contourf(data.lat, data.pfull, uv_conv_stat_after, extend = 'both', levels = np.arange(-5.,5.1,1.))
    ax4.contour(data.lat, data.pfull, u_ztav_after, colors='0.75', levels=np.arange(-60.,61.,10.), linewidths=2)
    ax4.contour(data.lat, data.pfull, psi_after, colors='k', levels=np.arange(-300.,301.,50.), linewidths=2)
    ax4.grid(True,linestyle=':')
    ax4.set_xlabel('Latitude')
    
    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
    #Colorbar
    cb1=f.colorbar(f2, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
    cb1.set_label('Eddy momentum flux convergence, $ms^{-1}day^{-1}$')
    
    plt.savefig(plot_dir+'upsi_eddies.pdf', format='pdf')
    plt.close()
    
    


upsi('ap_2', [121,361], [30,39])
upsi('full_qflux', [121,481], [18,39])
#upsi('flat_qflux', [121,481], [18,44], lonin = [60.,150.])


