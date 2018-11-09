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


def overturning_monthly(run, lonin=[-1.,361.]):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
        
    plot_dir = '/scratch/rg419/plots/overturning_monthly/' + run + '/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data.coords['month'] = (data.xofyear - 1) //6 + 1 
    data = data.groupby('month').mean(('xofyear'))
    #data['vcomp'] = data.vcomp.fillna(0.)
    #data['ucomp'] = data.ucomp.fillna(0.)
    # Create a VectorWind instance to handle the computation
    w = VectorWind(data.ucomp.sel(pfull=np.arange(50.,950.,50.)), data.vcomp.sel(pfull=np.arange(50.,950.,50.)))
    #w = VectorWind(data.ucomp, data.vcomp)
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    #print(vchi.pfull)
    #print(data.pfull)
    
    #data.vcomp.mean('lon')[0,:,:].plot.contourf(x='lat', y='pfull', yincrease=False, add_labels=False)
    #plt.figure(2)
    #vchi.mean('lon')[0,:,:].plot.contourf(x='lat', y='pfull', yincrease=False, add_labels=False)
    #plt.show()
    
    ds_chi = xr.Dataset({'vcomp': (vchi)},
                     coords={'month': ('month', vchi.month),
                             'pfull': ('pfull', vchi.pfull),
                               'lat': ('lat', vchi.lat),
                               'lon': ('lon', vchi.lon)})
    
    ds_psi = xr.Dataset({'vcomp': (vpsi)},
                     coords={'month': ('month', vchi.month),
                             'pfull': ('pfull', vchi.pfull),
                               'lat': ('lat', vchi.lat),
                               'lon': ('lon', vchi.lon)})
    
    def get_lons(lonin, data):
        if lonin[1]>lonin[0]:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
        else:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        return lons
    
    lons = get_lons(lonin,data)
        
    psi = mass_streamfunction(data, lons=lons, dp_in=50., use_v_locally=True)
    psi /= 1.e9
    
    psi_chi = mass_streamfunction(ds_chi, lons=lons, dp_in=50.)
    psi_chi /= 1.e9
    
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
        
        f1 = data.ucomp.sel(month=i+1).sel(lon=lons).mean('lon').plot.contourf(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-50.,50.1,5.), extend='both', add_labels=False, add_colorbar=False)
        
        m = mc.omega * mc.a**2. * np.cos(psi.lat*np.pi/180.)**2. + data.ucomp.sel(lon=lons).mean('lon') * mc.a * np.cos(psi.lat*np.pi/180.)
        m_levs = mc.omega * mc.a**2. * np.cos(np.arange(-60.,1.,5.)*np.pi/180.)**2.    
        m.sel(month=i+1).plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=m_levs, colors='0.7', add_labels=False)
        
        i=i+1
        ax.set_xlim(-35,35)
        ax.set_xticks(np.arange(-30,31,15))
        ax.grid(True,linestyle=':')
    
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Zonal wind speed, m/s')
    
    if lonin == [-1.,361.]:
        plt.savefig(plot_dir + 'psi_u_' + run + '.pdf', format='pdf')
    else:
        figname = plot_dir + 'psi_u_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
        plt.savefig(figname, format='pdf')
        
    plt.close()
    
    
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4)
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        psi_chi[:,i,:].plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(0.,601,100.), colors='k', add_labels=False)
        psi_chi[:,i,:].plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-600.,0.,100.), colors='k', linestyles='dashed', add_labels=False)
        
        f1 = uchi.sel(month=i+1).sel(lon=lons).mean('lon').plot.contourf(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-3.,3.1,0.5), extend='both', add_labels=False, add_colorbar=False)
        
        i=i+1
        ax.set_xlim(-35,35)
        ax.set_xticks(np.arange(-30,31,15))
        ax.grid(True,linestyle=':')
    
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Zonal wind speed, m/s')
    
    if lonin == [-1.,361.]:
        plt.savefig(plot_dir + 'psi_chi_u_' + run + '.pdf', format='pdf')
    else:
        figname = plot_dir + 'psi_chi_u_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
        plt.savefig(figname, format='pdf')
        
    plt.close()
    
    


overturning_monthly('half_shallow')
overturning_monthly('half_shallow', lonin=[350,10])
overturning_monthly('half_shallow', lonin=[80,100])
overturning_monthly('half_shallow', lonin=[170,190])
overturning_monthly('half_shallow', lonin=[260,280])

#overturning_monthly('half_shallow_5')
#overturning_monthly('half_shallow_5', lonin=[350,10])
#overturning_monthly('half_shallow_5', lonin=[80,100])
#overturning_monthly('half_shallow_5', lonin=[170,190])
#overturning_monthly('half_shallow_5', lonin=[260,280])

#overturning_monthly('half_shallow_10')
#overturning_monthly('half_shallow_10', lonin=[350,10])
#overturning_monthly('half_shallow_10', lonin=[80,100])
#overturning_monthly('half_shallow_10', lonin=[170,190])
#overturning_monthly('half_shallow_10', lonin=[260,280])

#overturning_monthly('q_shallow')
#overturning_monthly('q_shallow', lonin=[350,10])
#overturning_monthly('q_shallow', lonin=[35,45])
#overturning_monthly('q_shallow', lonin=[80,100])
#overturning_monthly('q_shallow', lonin=[215,235])

#overturning_monthly('3q_shallow')
#overturning_monthly('3q_shallow', lonin=[350,10])
#overturning_monthly('3q_shallow', lonin=[125,145])
#overturning_monthly('3q_shallow', lonin=[260,280])
#overturning_monthly('3q_shallow', lonin=[305,325])

