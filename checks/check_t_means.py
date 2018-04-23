#Load in u, v, t, w and check lat-lon pictures over a 5 year time average - any waves?

#check also monthly mean latitudinal phi structure

from time_av_xr import load_year_xr
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

#Initialise arrays to load into
run_fol = 'flat_10m_diags/np16'
rundata = load_year_xr(run_fol, 2, pinterp=True)
rundata.coords['month'] = (rundata.time // 30) - 11
mngrp = rundata.ucomp.groupby('month').mean(('time'))
u =     xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
v =     xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
w =     xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
phi =   xr.DataArray(np.zeros((12,17,64,128,5)), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])

for year in range(2,7):
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['month'] = (rundata.time // 30) + 1
    u[:,:,:,year-2]     = rundata.ucomp.mean(('time'))
    v[:,:,:,year-2]     = rundata.vcomp.mean(('time'))
    w[:,:,:,year-2]     = rundata.omega.mean(('time'))
    phi[:,:,:,:,year-2]   = rundata.height.groupby('month').mean(('time'))
    
u_av = u.mean(('year'))
v_av = v.mean(('year'))
w_av = w.mean(('year'))
phi_av = phi.mean(('year'))

plt.figure(0)
u_av[9,:,:].plot.contourf(x='lon', y='lat')

plt.figure(1)
v_av[9,:,:].plot.contourf(x='lon', y='lat')

plt.figure(2)
w_av[9,:,:].plot.contourf(x='lon', y='lat')

(phi_av[:,2,:,:]-phi_av[:,2,:,:].mean(('lon'))).plot.contourf(x='lon', y='lat',levels=np.arange(-180.,185.,10.),col='month', col_wrap=4)

plt.show()





