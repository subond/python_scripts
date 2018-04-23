"""
Function to plot meridional overturning, zonal wind speed, and eddy flux convergences

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means
from physics import mass_streamfunction, gradients as gr
import sh
from finite_difference import cfd
from pylab import rcParams

    
rcParams['figure.figsize'] = 12, 7
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/steady_state_runs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def uv_partitioning(data,lons):
        
    u_ztav = data.ucomp.sel(lon=lons).mean('lon') 
    v_ztav = data.vcomp.sel(lon=lons).mean('lon')
    
    uv_eddy = data.ucomp_vcomp.sel(lon=lons).mean('lon') - u_ztav * v_ztav 
                     
    uv_conv = gr.ddy(uv_eddy, uv=True) * -86400.

    return uv_conv, u_ztav
    
    
def u_psi(run, lonin=[-1.,361.]):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    psi = mass_streamfunction(data, a=6376.0e3, lons=lons, dp_in=50.)
    psi /= 1.e9    
    psi = psi.T

        
    uv_conv, u_ztav  = uv_partitioning(data, lons)
    
    ds = xr.Dataset({'psi': (['pfull', 'lat'], psi),
                     'uv_conv':  (['pfull', 'lat'], uv_conv),
                     'u_ztav':  (['pfull', 'lat'], u_ztav)},
                     coords={'pfull': ('pfull', data.pfull),
                               'lat': ('lat', data.lat)})
                               
    
                      
                               
    lwid=2

    ds.uv_conv.plot.contourf(x='lat', y='pfull', add_labels = False, levels = np.arange(-10.,10.1,1.), extend = 'both')

    plt.contour(ds.lat, ds.pfull, ds.u_ztav, colors='0.7', levels=np.arange(-60.,61.,10.), linewidths=lwid)
    psi.plot.contour(x='lat', y='pfull', colors='k', levels=np.arange(0.,301.,60.), linewidths=lwid)
    psi.plot.contour(x='lat', y='pfull', colors='k', linestyles='--', levels=np.arange(-300.,0.,60.), linewidths=lwid)
    plt.ylabel('Pressure, hPa')
    plt.xlabel('Latitude')
    plt.xlim(-60,60)
    plt.grid(True,linestyle=':')
    plt.gca().invert_yaxis()
    
    
    if lonin == [-1.,361.]:
        figname = 'psi_u_' + run + '.pdf'
    else:
        figname = 'psi_u_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
        
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    
    

u_psi('qflux_0_300', lonin=[50.,140.])
u_psi('qflux_0_300')
u_psi('qflux_5_300', lonin=[50.,140.])
u_psi('qflux_5_300')
u_psi('qflux_10_300', lonin=[50.,140.])
u_psi('qflux_10_300')
u_psi('qflux_15_300', lonin=[50.,140.])
u_psi('qflux_15_300')
u_psi('qflux_20_300', lonin=[50.,140.])
u_psi('qflux_20_300')
u_psi('qflux_25_300', lonin=[50.,140.])
u_psi('qflux_25_300')



