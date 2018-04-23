""" Investigate Emanuel 1995 critical subcloud entropy curvatures """

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from physics import gradients as gr



def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons
    
def subcloud_entropy_curve(data, lonin=[-1.,361.], rot_fac=1.):
    
    # declare constants
    cp = 287.04/2*7
    L = 2.500e6
    omega = 7.2921150e-5 * rot_fac
    a = 6376.0e3
    surf_trop_tdiff = 100.
    
    # pick longitudes to look at
    lons = pick_lons(data, lonin)
    
    # Evaluate equivalent potential temperature
    ept = ((data.temp + L/cp*data.sphum)*(1000./data.pfull)**(2./7.)).sel(pfull=850., lon=lons).mean('lon').drop('pfull')    
    data['ept'] = ept
    
    # moist entropy, s = cp ln(ept)
    s = cp*np.log(ept)
    
    # entropy curvature d/dy(cos(lat)^3/sin(lat) (Ts-Tt) ds/dy)
    dsdy = gr.ddy(s, vector=False)
    dsdy = dsdy * surf_trop_tdiff * np.cos(data.lat*np.pi/180.)**3. / np.sin(data.lat*np.pi/180.)
    data['s_curve'] = gr.ddy(dsdy, vector=False)
    
    # critical curvature = -4 * omega^2 * a^2 *cos(lat)^3 * sin(lat)
    crit_curve = -4. * omega**2 * a**2 * np.cos(data.lat*np.pi/180.)**3. * np.sin(data.lat*np.pi/180.)
    
    # Find max ept and location of max
    ept_max = [ept.max('lat'), ept.lat[ept.argmax('lat')]]
    
    # critical ept = ept_max * e^ (-chi * (cos(ept_max_lat)^2 - cos(lat)^2)^2 / cos(lat)^2)
    crit_ept_prefac = -1. * omega**2 * a**2/cp/surf_trop_tdiff
    
    n = len(data.xofyear.values)
    crit_ept_fn = np.zeros([n,64])
    crit_ept = np.zeros([n,64])
    
    for i in range(n):
        crit_ept_fn[i,:] = (np.cos(ept_max[1][i]*np.pi/180.)**2. - np.cos(data.lat*np.pi/180.)**2.)**2. / np.cos(data.lat*np.pi/180.)**2.
        crit_ept[i,:] = ept_max[0].values[i] * np.exp(crit_ept_prefac * crit_ept_fn[i,:])
        
    data['crit_ept'] = xr.DataArray(crit_ept, coords=[data.xofyear, data.lat], dims=['xofyear','lat'])
    
    # moist entropy, s = cp ln(ept)
    crit_s = cp*np.log(data.crit_ept)
    
    # entropy curvature d/dy(cos(lat)^3/sin(lat) (Ts-Tt) ds/dy)
    dcrit_sdy = gr.ddy(crit_s, vector=False)
    dcrit_sdy = dcrit_sdy * surf_trop_tdiff * np.cos(data.lat*np.pi/180.)**3. / np.sin(data.lat*np.pi/180.)
    data['crit_s_curve'] = gr.ddy(dcrit_sdy, vector=False)
    


def crit_curve_plot(run, lonin=[-1.,361.], rot_fac=1.):
    
    plot_dir = '/scratch/rg419/plots/subcloud_entropy/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    subcloud_entropy_curve(data, lonin=lonin, rot_fac=rot_fac)
    
    (data.s_curve[:,20:44] - data.crit_s_curve[:,20:44]).plot.contourf(x='xofyear', y='lat', levels=np.arange(-5.e-7, 5.1e-7, 5e-8))
    try:
        (data.precipitation*86400.).mean('lon').plot.contour(x='xofyear', y='lat', colors='k', levels = np.arange(2.,15.,2.))
    except:
        data['precipitation'] = (data.convection_rain + data.condensation_rain)*86400.
        data.precipitation.mean('lon').plot.contour(x='xofyear', y='lat', colors='k', levels = np.arange(2.,15.,2.))
    
    plt.savefig(plot_dir + 'crit_curve_' + run + '.pdf', format='pdf')
    plt.close()
    

if __name__ == "__main__":

    crit_curve_plot('zs_sst')    
    #crit_curve_plot('sn_1.000')
    #crit_curve_plot('sn_2.000')
    #crit_curve_plot('sn_0.500')
    #crit_curve_plot('rt_0.500', rot_fac=0.5)
    #crit_curve_plot('rt_2.000', rot_fac=2.)
    #crit_curve_plot('ap_2')
    #crit_curve_plot('ap_20')
    #crit_curve_plot('full_qflux', lonin=[60.,150.])
    
    
    
    #ept[:,20:44].plot.contourf(x='xofyear', y='lat', levels=np.arange(290.,381.,5.))
    #data.omega.mean('lon').sel(pfull=500)[:,20:44].plot.contour(x='xofyear', y='lat', colors='k')


    #plt.show()