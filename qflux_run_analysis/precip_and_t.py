# Plot the surface temperature and precipitation through the year for the monsoon region

from data_handling import time_means, month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import gc

plt.rc('text', usetex=True)
font = {'family' : 'sans-serif','sans-serif':['Helvetica'],
        'weight' : 'bold',
        'size'   : 18}

plt.rc('font', **font)

data = time_means('aquaplanet_2m', [121,481], filename='atmos_daily', timeav='pentad')


def plot_precip(data, lons=range(0,128)):
    print 'precip'
    data['totp'] = (('xofyear','lat','lon'), (data.convection_rain + data.condensation_rain)*86400.)

    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    ax=data.totp[:,:,lons].mean('lon').plot.contourf(x='xofyear', y='lat',levels=np.arange(6.,25.,3.), add_label = False, add_colorbar=False, extend='min')
    cb1=plt.colorbar(ax)
    cb1.set_label('Precip, mm/day')
    cs = data.t_surf[:,:,lons].mean('lon').plot.contour(x='xofyear', y='lat',levels=range(250,321,10), add_label = False, colors='w', add_colorbar=False)
    plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
    plt.ylim((-45,45))
    plt.xlim((1,73))
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.ylabel('Latitude')
    #plt.title('Zonal mean precipitation evolution, mm/day')
    plt.tight_layout()


def plot_u(data, lons=range(0,128)):
    print 'u'
    data['ucomp_850'] = (('xofyear','lat','lon'), data.ucomp[:,2,:,:].values)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    ax = data.ucomp_850[:,:,lons].mean('lon').plot.contourf(x='xofyear', y='lat',levels=np.arange(-20.,21.,4.), add_label = False, add_colorbar=False, extend='neither')
    cb1=plt.colorbar(ax)
    cb1.set_label('Zonal wind speed, m/s (850 hPa)')
    cs = data.t_surf[:,:,lons].mean('lon').plot.contour(x='xofyear', y='lat',levels=range(250,321,10), add_label = False, colors='w', add_colorbar=False)
    plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
    plt.ylim((-45,45))
    plt.xlim((1,73))
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.ylabel('Latitude')
    #plt.title('Zonal mean precipitation evolution, mm/day')
    plt.tight_layout()

schneider = [i for i in range(len(data.lon)) if data.lon[i] >= 70. and data.lon[i] < 100.]
africa = [i for i in range(len(data.lon)) if data.lon[i] >= 330. or data.lon[i] < 60.]
asia = [i for i in range(len(data.lon)) if data.lon[i] >= 60. and data.lon[i] < 150.]
pacific = [i for i in range(len(data.lon)) if data.lon[i] >= 150. and data.lon[i] < 240.]
america = [i for i in range(len(data.lon)) if data.lon[i] >= 240. and data.lon[i] < 330.]

plot_precip(data)
plt.savefig('/scratch/rg419/plots/qflux_run_analysis/precip_ap2.png')
plt.clf()

#plot_u(data,schneider)
#plt.savefig('/scratch/rg419/plots/qflux_run_analysis/u_bs.png')
#plt.clf()
#plot_precip(data,schneider)
#plt.savefig('/scratch/rg419/plots/qflux_run_analysis/precip_bs.png')
#plt.clf()


