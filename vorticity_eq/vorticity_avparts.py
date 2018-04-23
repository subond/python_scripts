# Evaluate 'leading order' terms in vorticity equation: mean horizontal vel dotted with mean gradient of abs vort, and mean abs vort times mean vel divergence

from physics import mass_streamfunction
from data_handling import time_means
from finite_difference import cfd
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np


def vort_av(run, months, filename='atmos_pentad', timeav='pentad', period_fac=1.):
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    abs_vort = data.vor + f
    
    dvortdx = xr.DataArray( cfd(abs_vort.values,data.lon*np.pi/180.,3),   [('xofyear', data.xofyear ), ('pfull', data.pfull ), ('lat', data.lat), ('lon', data.lon )])
    dvortdx = dvortdx/coslat/a
    dvortdy = xr.DataArray( cfd(abs_vort.values,data.lat*np.pi/180.,2),   [('xofyear', data.xofyear ), ('pfull', data.pfull ), ('lat', data.lat), ('lon', data.lon )])
    dvortdy = dvortdy/a
    dvortdp = xr.DataArray( cfd(abs_vort.values,data.pfull*100.,1),   [('xofyear', data.xofyear ), ('pfull', data.pfull ), ('lat', data.lat), ('lon', data.lon )])
    
    vort_adv = data.ucomp*dvortdx + data.vcomp*dvortdy
    vort_vert = data.omega*dvortdp
    vort_div = abs_vort*data.div
    
    precip = (data.convection_rain + data.condensation_rain)*86400.
    
    plt.figure(1)
    vort_adv[:,15,:,:].mean('lon').plot.contourf(x='xofyear', y='lat',levels = np.arange(-4.5e-10, 4.6e-10, 0.5e-10))
    plt.figure(2)
    vort_div[:,15,:,:].mean('lon').plot.contourf(x='xofyear', y='lat',levels = np.arange(-4.5e-10, 4.6e-10, 0.5e-10))
    plt.figure(3)
    vort_vert[:,15,:,:].mean('lon').plot.contourf(x='xofyear', y='lat')
    plt.figure(4)
    abs_vort[:,15,:,:].mean('lon').plot.contourf(x='xofyear', y='lat',levels = np.arange(-2.4e-4, 2.5e-4, 0.2e-4))
    precip.mean('lon').plot.contour(x='xofyear', y='lat',levels=np.arange(6.,31.,6.), add_label = False, add_colorbar=False,colors='k')
    plt.show()
    
    
    
vort_av('ap_1_partraw', [193,313])