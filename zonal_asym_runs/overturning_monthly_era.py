''' 
03/08/2018 Plot local mass streamfunction from climatology month by month
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import mass_streamfunction
from data_handling_updates import model_constants as mc
from windspharm.xarray import VectorWind


def overturning_monthly(lonin=[-1.,361.]):
    
    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_v_clim_alllevs.nc' )
    data_u = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_u_clim_alllevs.nc' )
        
    plot_dir = '/scratch/rg419/plots/overturning_monthly/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    def get_lons(lonin, data):
        if lonin[1]>lonin[0]:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
        else:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        return lons
    
    lons = get_lons(lonin,data)
    
    month_lengths = [31,28,31,30,31,30,31,31,30,31,30,31]
    month = []
    for i in range(1,13):
        month = month + [i]*month_lengths[i-1]
        
    data.coords['month'] = data.day_of_yr*0. + np.array(month)
    data = data.groupby('month').mean('day_of_yr')   
    
    data_u.coords['month'] = data_u.day_of_yr*0. + np.array(month)
    data_u = data_u.groupby('month').mean('day_of_yr')
    
    data = data.rename(name_dict={'v': 'vcomp'})
    data_u = data_u.rename(name_dict={'u': 'ucomp'})
    
    # Create a VectorWind instance to handle the computation
    w = VectorWind(data_u.ucomp.sel(pfull=np.arange(50.,950.,50.)), data.vcomp.sel(pfull=np.arange(50.,950.,50.)))
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    ds_chi = xr.Dataset({'vcomp': (vchi)},
                     coords={'month': ('month', vchi.month),
                             'pfull': ('pfull', vchi.pfull),
                               'lat': ('lat', vchi.lat),
                               'lon': ('lon', vchi.lon)})
    
    psi = mass_streamfunction(data, lons=lons, dp_in=-50., intdown=False)
    psi /= 1.e9
    
    psi_chi = mass_streamfunction(ds_chi, lons=lons, dp_in=-50., intdown=False)
    psi_chi /= 1.e9
    
   # print(psi)
    # Set figure parameters
    rcParams['figure.figsize'] = 10, 7
    rcParams['font.size'] = 14
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4)
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        psi[:,i,:].plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(0.,601,100.), colors='k', add_labels=False)
        psi[:,i,:].plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-600.,0.,100.), colors='k', linestyles='dashed', add_labels=False)
        
        f1 = data_u.ucomp.sel(month=i+1).sel(lon=lons).mean('lon').plot.contourf(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-50.,50.1,5.), extend='both', add_labels=False, add_colorbar=False)
        
        m = mc.omega * mc.a**2. * np.cos(psi.lat*np.pi/180.)**2. + data_u.ucomp.sel(lon=lons).mean('lon') * mc.a * np.cos(psi.lat*np.pi/180.)
        m_levs = mc.omega * mc.a**2. * np.cos(np.arange(-60.,1.,5.)*np.pi/180.)**2.    
        m.sel(month=i+1).plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=m_levs, colors='0.7', add_labels=False)
        
        i=i+1
        ax.set_xlim(-35,35)
        ax.set_xticks(np.arange(-30,31,15))
        ax.grid(True,linestyle=':')
    
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    #cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    #cb1.set_label('Zonal wind speed, m/s')
    
    if lonin == [-1.,361.]:
        plt.savefig(plot_dir + 'psi_u_era.pdf', format='pdf')
    else:
        figname = plot_dir + 'psi_u_era_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
        plt.savefig(figname, format='pdf')
        
    plt.close()
    
    
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4)
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        psi_chi[:,i,:].plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(0.,601,100.), colors='k', add_labels=False)
        psi_chi[:,i,:].plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-600.,0.,100.), colors='k', linestyles='dashed', add_labels=False)
                
        i=i+1
        ax.set_xlim(-35,35)
        ax.set_xticks(np.arange(-30,31,15))
        ax.grid(True,linestyle=':')
    
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Zonal wind speed, m/s')
    
    if lonin == [-1.,361.]:
        plt.savefig(plot_dir + 'psi_chi_era.pdf', format='pdf')
    else:
        figname = plot_dir + 'psi_chi_era_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
        plt.savefig(figname, format='pdf')
        
    plt.close()


overturning_monthly()
overturning_monthly(lonin=[100,150])
overturning_monthly(lonin=[70,100])
overturning_monthly(lonin=[15,60])
