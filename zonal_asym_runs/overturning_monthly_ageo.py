''' 
03/08/2018 Plot local mass streamfunction from climatology month by month
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from hadley_cell import mass_streamfunction
from data_handling_updates import model_constants as mc, gradients as gr
from windspharm.xarray import VectorWind


def overturning_monthly(run, lonin=[-1.,361.]):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
        
    plot_dir = '/scratch/rg419/plots/overturning_monthly_ageo/' + run + '/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data.coords['month'] = (data.xofyear - 1) //6 + 1 
    data = data.groupby('month').mean(('xofyear'))
    
    #Coriolis
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    dphidx = gr.ddx(data.height)
    dphidx = 9.8 * dphidx
    v_ageo = data.vcomp - dphidx/f
    v_geo = dphidx/f
    
    dphidy = gr.ddy(data.height, vector=False)
    dphidy = 9.8 * dphidy
    u_ageo = data.ucomp - dphidy/f
    u_geo = dphidx/f
    
    # Create a VectorWind instance to handle the computation
    #w = VectorWind(data.ucomp.sel(pfull=np.arange(50.,950.,50.)), data.vcomp.sel(pfull=np.arange(50.,950.,50.)))
    w = VectorWind(u_ageo.sel(pfull=np.arange(50.,950.,50.)), v_ageo.sel(pfull=np.arange(50.,950.,50.)))
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    
    ds_chi = xr.Dataset({'vcomp': (vchi)},
                     coords={'month': ('month', vchi.month),
                             'pfull': ('pfull', vchi.pfull),
                               'lat': ('lat', vchi.lat),
                               'lon': ('lon', vchi.lon)})
    
    w = VectorWind(u_geo.sel(pfull=np.arange(50.,950.,50.)), v_geo.sel(pfull=np.arange(50.,950.,50.)))
    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    
    ds_chi_geo = xr.Dataset({'vcomp': (vchi)},
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
    
    
    psi = mass_streamfunction(data, lons=lons, dp_in=50.)
    psi /= 1.e9
    
    psi_chi = mass_streamfunction(ds_chi, lons=lons, dp_in=-50., intdown=False)
    psi_chi /= 1.e9
    
    psi_chi_geo = mass_streamfunction(ds_chi_geo, lons=lons, dp_in=-50., intdown=False)
    psi_chi_geo /= 1.e9
    
    # Set figure parameters
    rcParams['figure.figsize'] = 10, 7
    rcParams['font.size'] = 14
    
    
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
    
    if lonin == [-1.,361.]:
        plt.savefig(plot_dir + 'psi_chi_' + run + '.pdf', format='pdf')
    else:
        figname = plot_dir + 'psi_chi_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
        plt.savefig(figname, format='pdf')
        
    plt.close()
    
    
    # Start figure with 12 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4)
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
    
    i=0
    for ax in axes:
        psi_chi_geo[:,i,:].plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(0.,601,100.), colors='k', add_labels=False)
        psi_chi_geo[:,i,:].plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-600.,0.,100.), colors='k', linestyles='dashed', add_labels=False)
        
        i=i+1
        ax.set_xlim(-35,35)
        ax.set_xticks(np.arange(-30,31,15))
        ax.grid(True,linestyle=':')
    
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    if lonin == [-1.,361.]:
        plt.savefig(plot_dir + 'psi_chi_geo_' + run + '.pdf', format='pdf')
    else:
        figname = plot_dir + 'psi_chi_geo_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
        plt.savefig(figname, format='pdf')
        
    plt.close()
    
    


#overturning_monthly('half_shallow')
overturning_monthly('half_shallow', lonin=[350,10])
overturning_monthly('half_shallow', lonin=[80,100])
overturning_monthly('half_shallow', lonin=[170,190])
overturning_monthly('half_shallow', lonin=[260,280])

overturning_monthly('half_shallow_5')
overturning_monthly('half_shallow_5', lonin=[350,10])
overturning_monthly('half_shallow_5', lonin=[80,100])
overturning_monthly('half_shallow_5', lonin=[170,190])
overturning_monthly('half_shallow_5', lonin=[260,280])

overturning_monthly('half_shallow_10')
overturning_monthly('half_shallow_10', lonin=[350,10])
overturning_monthly('half_shallow_10', lonin=[80,100])
overturning_monthly('half_shallow_10', lonin=[170,190])
overturning_monthly('half_shallow_10', lonin=[260,280])

#overturning_monthly('q_shallow')
#overturning_monthly('q_shallow', lonin=[350,10])
#overturning_monthly('q_shallow', lonin=[35,45])
#overturning_monthly('q_shallow', lonin=[80,100])
#overturning_monthly('q_shallow', lonin=[215,235])

overturning_monthly('3q_shallow')
overturning_monthly('3q_shallow', lonin=[350,10])
overturning_monthly('3q_shallow', lonin=[125,145])
overturning_monthly('3q_shallow', lonin=[260,280])
overturning_monthly('3q_shallow', lonin=[305,325])

