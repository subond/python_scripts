from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from plotting_fns import plot_var_ll, plot_quiver_ll, plot_var_zav
from time_av_fns import multi_month_load, monthly_mean, intermonthly_mean

#set run name
run_fol = 'land_test'

#open files

nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/' + run_fol + '/np8/run1/atmos_daily.nc'
fh = Dataset(nc_file, mode='r')
land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/' + run_fol + '/input/land.nc'
lfh = Dataset(land_file, 'r', format='NETCDF3_CLASSIC')
land_mask = lfh.variables['land_mask'][:]
lfh.close()

#read in variables
lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
p = fh.variables['pfull'][:]
times = fh.variables['time'][:]
fh.close()

#temp = fh.variables['temp'][:]
#t_surf = fh.variables['t_surf'][:]
#slp = fh.variables['slp'][:]
#u = fh.variables['ucomp'][:]
#v = fh.variables['vcomp'][:]
#w = fh.variables['omega'][:]
#z = fh.variables['height'][:]
#rh = fh.variables['rh'][:]
#cnvp = fh.variables['convection_rain'][:]
#cndp = fh.variables['condensation_rain'][:]


u = intermonthly_mean('ucomp',run_fol, [6,7,8])
v = intermonthly_mean('vcomp',run_fol, [6,7,8])
temp = intermonthly_mean('temp',run_fol, [6,7,8])
tsurf = intermonthly_mean('t_surf',run_fol, [6,7,8],twodim=True)
cnvp = intermonthly_mean('convection_rain',run_fol, [6,7,8],twodim=True)
cndp = intermonthly_mean('condensation_rain',run_fol, [6,7,8],twodim=True)

plot_var_ll(lons,lats,tsurf,land_mask=land_mask, plotno=1,pltitle = 'Surface Temperature', cblabel = 'K')

plot_var_ll(lons,lats,cnvp*86400,land_mask=land_mask, plotno=5,pltitle = 'Convective Precip', cblabel = 'mm/day')

plot_var_ll(lons,lats,cndp*86400,land_mask=land_mask, plotno=6,pltitle = 'Large-scale Precip', cblabel = 'mm/day')

plot_var_ll(lons,lats,(cndp + cnvp)*86400,land_mask=land_mask, plotno=7,pltitle = 'Total Precip', cblabel = 'mm/day')

plot_quiver_ll(lons,lats,u[26,:,:],v[26,:,:],land_mask=land_mask, plotno=4,pltitle = 'Horizontal windspeed, m/s')

plot_var_zav(lats,p,np.mean(u,2),v=range(-25,55,5),plotno=2)
plot_var_zav(lats,p,np.mean(temp,2),v=range(180,320,10), plotno=3) 

plt.show()
