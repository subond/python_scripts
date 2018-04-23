#load up diags from the run including the bucket model:
#1) check if it had land in it...
#2) check initial state
#3) check contributions to depth are consistent with depth at different times
#4) look at behaviour over the seasonal cycle: change in bucket depth and latent heat release

import sys
sys.path.append('/scratch/rg419/workdir_moist/python_bin')
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from time_av_xr import load_year_xr

#set run name
run_fol = 'bucket_test2/np8'
inp_fol = 'bucket_test2'

land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/' + inp_fol + '/input/land.nc'
land = xr.open_dataset( land_file)


rundata = load_year_xr(run_fol, [1])

rundata.coords['month'] = (rundata.time // 30) + 1

#take time and zonal mean
bd_tav = rundata.bucket_depth.groupby('month').mean(('time'))
dbdlh_tav = rundata.bucket_depth_lh.groupby('month').mean(('time'))
dbdcd_tav = rundata.bucket_depth_cond.groupby('month').mean(('time'))
dbdcv_tav = rundata.bucket_depth_conv.groupby('month').mean(('time'))

raincv_tav = rundata.convection_rain.groupby('month').mean(('time'))
raincd_tav = rundata.condensation_rain.groupby('month').mean(('time'))
lhe_tav = rundata.flux_lhe.groupby('month').mean(('time'))
totrain_tav = (raincv_tav + raincd_tav)*86400.

rundata_month = rundata.groupby('month').mean(('time'))

depth_change = rundata.bucket_depth_cond + rundata.bucket_depth_conv - rundata.bucket_depth_lh
depth_change_bucket = rundata.bucket_depth[0:30,:,:].diff(('time'))

isotach = np.sqrt(rundata_month.ucomp[:,38,:,:]**2 + rundata_month.vcomp[:,38,:,:]**2)


rundata_month.bucket_depth.plot.contourf(x='lon', y='lat',levels=np.arange(0,2.,0.1), col='month', col_wrap=4, add_label = False)
plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/bucket_depth_land.png')

rundata_month.bucket_depth.plot.contourf(x='lon', y='lat',levels=np.arange(0,24.,1.), col='month', col_wrap=4, add_label = False)
plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/bucket_depth_sea.png')

totrain_tav.plot.contourf(x='lon', y='lat',levels=np.arange(0,20.,0.5), col='month', col_wrap=4, add_label = False)
plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/tot_precip.png')

rundata_month.flux_lhe.plot.contourf(x='lon', y='lat',levels=np.arange(0,280.,10.), col='month', col_wrap=4, add_label = False)
plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/evap.png')

rundata_month.omega[:,37,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.2,0.21,0.01), col='month', col_wrap=4, add_label = False)
plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/omega.png')

rundata_month.t_surf.plot.contourf(x='lon', y='lat',levels=np.arange(230,321,5), col='month', col_wrap=4, add_label = False)
plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/t_surf.png')

isotach.plot.contourf(x='lon', y='lat',levels=np.arange(0.,32.,1.), col='month', col_wrap=4, add_label = False)

plt.show()









#dbdlh_tav.plot.contourf(x='lon', y='lat', col='month', col_wrap=4, add_label = False)
#dbdcd_tav.plot.contourf(x='lon', y='lat', col='month', col_wrap=4, add_label = False)
#dbdcv_tav.plot.contourf(x='lon', y='lat', col='month', col_wrap=4, add_label = False)
#raincv_tav.plot.contourf(x='lon', y='lat', col='month', col_wrap=4, add_label = False)
#raincd_tav.plot.contourf(x='lon', y='lat', col='month', col_wrap=4, add_label = False)
#rundata.bucket_depth[0:30,:,:].plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#rundata.bucket_depth_lh[0:30,:,:].plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#rundata.bucket_depth_cond[0:30,:,:].plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#rundata.bucket_depth_conv[0:30,:,:].plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#depth_change[0:30,:,:].plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#depth_change_bucket.plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#plt.show()

#rundata.bucket_depth[30:60,:,:].plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#rundata.bucket_depth_lh[30:60,:,:].plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#rundata.bucket_depth_cond[30:60,:,:].plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#rundata.bucket_depth_conv[30:60,:,:].plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#depth_change[30:60,:,:].plot.contourf(x='lon', y='lat', col='time', col_wrap=6, add_label = False)
#plt.show()
