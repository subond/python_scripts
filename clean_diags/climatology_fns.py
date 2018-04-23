# Produce plots of the zonal mean terms in the momentum budget as a function of time for the aquaplanet simulations, showing the transition and change in balance of terms

from data_handling import time_means, season_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from physics import mass_streamfunction


def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons
    
def rain_and_sst(run, months, filename='plev_pentad', period_fac=1., plot_land = True):
    rcParams['figure.figsize'] = 20, 10
    
    plot_dir = '/scratch/rg419/plots/climatology/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    sn_dic = season_dic(1)
    
    data = time_means(run,months, filename=filename, timeav='season', period_fac=period_fac)
    data['rain'] = (('xofyear','lat','lon'), (data.convection_rain + data.condensation_rain)*86400.)

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    data['land'] = (('lat','lon'), land.land_mask)
 
    for i in range(0,4):
        print i
        ax = data.rain[i,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(0.,31.,3.), add_labels = False, extend='max', add_colorbar=False)
        if plot_land:
            data.land.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), colors='k', add_colorbar=False, add_labels=False)
        cs = data.t_surf[i,:,:].plot.contour(x='lon', y='lat', levels = np.arange(200., 350., 10.), colors='w', add_labels = False, add_colorbar=False)
        plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
        cb1=plt.colorbar(ax)
        cb1.set_label('Precipitation, mm/day')
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.tight_layout()  
        plt.savefig(plot_dir+'rain_and_sst_'+ sn_dic[i] + '.png')
        plt.close()    


def wind_and_isotachs(run, months, filename='plev_pentad', period_fac=1., plot_land = True):
    rcParams['figure.figsize'] = 20, 10
    
    plot_dir = '/scratch/rg419/plots/climatology/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    sn_dic = season_dic(1)
    
    data = time_means(run,months, filename=filename, timeav='season', period_fac=period_fac)
    data['isotach'] = (('xofyear','pfull','lat','lon'), np.sqrt(data.ucomp**2. + data.vcomp**2.))

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    data['land'] = (('lat','lon'), land.land_mask)
    
    
    for i in range(0,4):
        print i
        ax = data.isotach[i,3,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(0.,20.,2.), add_labels = False, add_colorbar=False, extend='max')
        if plot_land:
            data.land.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), colors='k', add_colorbar=False, add_labels=False)
        plt.quiver(data.lon[::3], data.lat[::1], data.ucomp[i,4,::1,::3], data.vcomp[i,4,::1,::3], headlength=3, headwidth=2)#, scale=500.,headwidth=5)
        cb1=plt.colorbar(ax)
        cb1.set_label('850 hPa wind speed, m/s')
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.tight_layout()  
        plt.savefig(plot_dir+'wind_and_isotachs_850_'+ sn_dic[i] + '.png')
        plt.close()    

    for i in range(0,4):
        print i
        ax = data.isotach[i,16,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(0.,60.,5.), add_labels = False, add_colorbar=False, extend='max')
        if plot_land:
            data.land.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), colors='k', add_colorbar=False, add_labels=False)
        plt.quiver(data.lon[::3], data.lat[::1], data.ucomp[i,16,::1,::3], data.vcomp[i,16,::1,::3], headlength=3, headwidth=2)#, scale=500.,headwidth=5)
        cb1=plt.colorbar(ax)
        cb1.set_label('200 hPa wind speed, m/s')
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.tight_layout()  
        plt.savefig(plot_dir+'wind_and_isotachs_200_'+ sn_dic[i] + '.png')
        plt.close()
        
            
        
def surface_fluxes(run, months, filename='plev_pentad', period_fac=1., plot_land = True):
    rcParams['figure.figsize'] = 20, 10
    stefan = 5.6734e-8
    
    plot_dir = '/scratch/rg419/plots/climatology/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    sn_dic = season_dic(1)
    
    data = time_means(run,months, filename=filename, timeav='season', period_fac=period_fac)
    data['surf_flux'] = (('xofyear','lat','lon'), stefan*data.t_surf**4. - data.flux_lw - data.flux_lw + data.flux_t + data.flux_lhe)

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    data['land'] = (('lat','lon'), land.land_mask)
 
    for i in range(0,4):
        print i
        ax = data.surf_flux[i,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-700.,700.,50.), add_labels = False, extend='both', add_colorbar=False)
        if plot_land:
            data.land.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), colors='k', add_colorbar=False, add_labels=False)
        cb1=plt.colorbar(ax)
        cb1.set_label('Surface flux, W/m2')
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.tight_layout()  
        plt.savefig(plot_dir+'surface_fluxes_'+ sn_dic[i] + '.png')
        plt.close() 


def diabatic_heating(run, months, filename='plev_pentad', period_fac=1., plot_land = True):
    rcParams['figure.figsize'] = 20, 10
    stefan = 5.6734e-8
    
    plot_dir = '/scratch/rg419/plots/climatology/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    sn_dic = season_dic(1)
    
    data = time_means(run,months, filename=filename, timeav='season', period_fac=period_fac)
    data['heat_rate'] = (('xofyear','pfull','lat','lon'), (data.tdt_rad + data.dt_tg_condensation + data.dt_tg_convection + data.dt_tg_diffusion)*86400.)

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc'
    land = xr.open_dataset( land_file)
    data['land'] = (('lat','lon'), land.land_mask)
 
    for i in range(0,4):
        print i
        ax = data.heat_rate[i,4,:,:].plot.contourf(x='lon', y='lat', levels=np.arange(-10.,10.,1.), add_labels = False, extend='both', add_colorbar=False)
        if plot_land:
            data.land.plot.contour(x='lon', y='lat', levels=np.arange(0.,2.,1.), colors='k', add_colorbar=False, add_labels=False)
        cb1=plt.colorbar(ax)
        cb1.set_label('Diabatic heating, K/day')
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.tight_layout()  
        plt.savefig(plot_dir+'diabatic_heating_'+ sn_dic[i] + '.png')
        plt.close() 
        


def u_and_t_zmean(run, months, filename='plev_pentad', period_fac=1., lonin=[-1.,361.]):
    rcParams['figure.figsize'] = 12, 10
    
    plot_dir = '/scratch/rg419/plots/climatology/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    sn_dic = season_dic(1)
    
    data = time_means(run,months, filename=filename, timeav='season', period_fac=period_fac)
    lons = pick_lons(data, lonin)
    
    uwnd = data.ucomp[:,:,:,lons].mean('lon')
    temp = data.temp[:,:,:,lons].mean('lon')
    
    for i in range(0,4):
        print i
        ax = uwnd[i,:,:].plot.contourf(x='lat', y='pfull', levels=np.arange(-60.,60.,5.), add_label = False, add_colorbar=False, yincrease=False, extend='both')
        cs = temp[i,:,:].plot.contour(x='lat', y='pfull', levels=np.arange(200.,300.,10.), colors='k', add_label = False, add_colorbar=False, yincrease=False, extend='both')
        plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
        cb1=plt.colorbar(ax)
        cb1.set_label('Zonal wind speed, m/s')
        plt.ylabel('Pressure, hPa')
        plt.xlabel('Latitude')
        plt.tight_layout()  
        plt.savefig(plot_dir+'u_and_t_zmean_'+ sn_dic[i] + '.png')
        plt.close()    



def u_and_psi_zmean(run, months, filename='atmos_pentad', period_fac=1., lonin=[-1.,361.]):
    rcParams['figure.figsize'] = 12, 10
    
    plot_dir = '/scratch/rg419/plots/climatology/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    sn_dic = season_dic(1)
    
    data = time_means(run,months, filename=filename, timeav='season', period_fac=period_fac)
    lons = pick_lons(data, lonin)
    
    uwnd = data.ucomp[:,:,:,lons].mean('lon')
    
    psi = mass_streamfunction(data, a=6376.0e3, lons=lons, dp_in=50.)
    psi /= 1.e9
        
    for i in range(0,4):
        print i
        ax = uwnd[i,:,:].plot.contourf(x='lat', y='pfull', levels=np.arange(-60.,60.,5.), add_label = False, add_colorbar=False, yincrease=False, extend='both')
        cs = psi[:,i,:].plot.contour(x='lat', y='pfull', levels=np.arange(-350.,350.,50.), colors='k', add_label = False, add_colorbar=False, yincrease=False, extend='both')
        plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
        cb1=plt.colorbar(ax)
        cb1.set_label('Zonal wind speed, m/s')
        plt.ylabel('Pressure, hPa')
        plt.xlabel('Latitude')
        plt.tight_layout()  
        plt.savefig(plot_dir+'u_and_psi_zmean_alllons_plev'+ sn_dic[i] + '.png')
        plt.close()    
        


#surface_fluxes('ap_1_rd', [121,481], plot_land = False)
#diabatic_heating('ap_1_rd', [121,481], plot_land = False)
#rain_and_sst('ap_1_rd', [121,481], plot_land = False)
#wind_and_isotachs('ap_1_rd', [121,481], plot_land = False)
#u_and_t_zmean('ap_1_rd', [121,481])
#u_and_psi_zmean('ap_1_rd', [121,481])

#surface_fluxes('ap_20', [121,481], plot_land = False)
#diabatic_heating('ap_20', [121,481], plot_land = False)
#rain_and_sst('ap_20', [121,481], plot_land = False)
#wind_and_isotachs('ap_20', [121,481], plot_land = False)
#u_and_t_zmean('ap_20', [121,481])
#u_and_psi_zmean('ap_20', [121,481])

#surface_fluxes('full_qflux', [121,481])
#diabatic_heating('full_qflux', [121,481])
#rain_and_sst('full_qflux', [121,481])
#wind_and_isotachs('full_qflux', [121,481])
#u_and_t_zmean('full_qflux', [121,481], lonin=[60.,150.])
u_and_psi_zmean('full_qflux', [121,481],filename='plev_pentad')# lonin=[60.,150.])

#surface_fluxes('flat_qflux', [121,397])
#diabatic_heating('flat_qflux', [121,397])
#rain_and_sst('flat_qflux', [121,397])
#wind_and_isotachs('flat_qflux', [121,397])
#u_and_t_zmean('flat_qflux', [121,397], lonin=[60.,150.])
u_and_psi_zmean('flat_qflux', [121,481],filename='plev_pentad')# lonin=[60.,150.])
