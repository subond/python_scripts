# Plot the surface temperature and precipitation through the year for the monsoon region

from data_handling import time_means, month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh

font = {'size'   : 18}
plt.rc('font', **font)

def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons
    
def p_and_t(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1., lonin=[-1.,361.]):
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    totp = (data.convection_rain + data.condensation_rain)*86400.
    lons = pick_lons(data, lonin)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]    
    
#    ax=totp[:,:,lons].mean('lon').plot.contourf(x='xofyear', y='lat',levels=np.arange(6.,31.,3.), add_label = False, add_colorbar=False)
    ax=totp[:,:,lons].mean('lon').plot.contourf(x='xofyear', y='lat',levels=np.arange(2.,15.,2.), add_label = False, add_colorbar=False, extend='both')
    totp[:,:,lons].mean('lon').plot.contour(x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_label = False, add_colorbar=False, colors='k')
    plt.grid(True,linestyle=':')
    cb1=plt.colorbar(ax)
    cb1.set_label('Precip, mm/day')
    cs = data.t_surf[:,:,lons].mean('lon').plot.contour(x='xofyear', y='lat',levels=range(250,321,10), add_label = False, colors='w', add_colorbar=False)
    plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
    plt.ylim((-45,45))
    plt.xlim((1,73))
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.ylabel('Latitude')
    plt.tight_layout()  
    plt.savefig(plot_dir+ 'p_and_t.png')
    plt.close()


#p_and_t('full_qflux',[121,481],lonin=[60.,150.])
#p_and_t('flat_qflux',[121,481],lonin=[60.,150.])
p_and_t('ap_2',[121,481])
p_and_t('am_qflux',[121,481],lonin=[60.,150.])
p_and_t('ap_20_qflux',[121,481],lonin=[60.,150.])
