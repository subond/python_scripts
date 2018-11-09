'''5/09/2018 Use JRA-55 and GPCP monthly data to plot 1979-2016 wind and precip difference between MJJAS and NDJFM'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams

plot_dir = '/scratch/rg419/plots/monsoon_review_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
# Load in data
data_u = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/ucomp_monthly/atmos_monthly_together.nc')
data_v = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_monthly/atmos_monthly_together.nc')
name_temp = '/scratch/rg419/obs_and_reanalysis/GPCP_monthly/gpcp_cdr_v23rB1_y%04d_m%02d.nc'
names = [name_temp % (y, m) for y in range(1979,2017) for m in range(1,13)]
data_p = xr.open_mfdataset(names)
land_mask='/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
land = xr.open_dataset(land_mask)

print(data_u.lon)
print(data_p.longitude)
print(land.longitude)

# Make climatologies
u_clim = data_u.groupby('time.month').mean('time')
v_clim = data_v.groupby('time.month').mean('time')
p_clim = data_p.groupby('time.month').mean('time')
print('means taken')

# Calculate MJJAS and NDJFM difference
u_diff = u_clim.sel(month=[5,6,7,8,9]).mean('month') - u_clim.sel(month=[11,12,1,2,3]).mean('month')
v_diff = v_clim.sel(month=[5,6,7,8,9]).mean('month') - v_clim.sel(month=[11,12,1,2,3]).mean('month')
p_diff = p_clim.sel(month=[5,6,7,8,9]).mean('month') - p_clim.sel(month=[11,12,1,2,3]).mean('month')

# Define monsoon region using precip
monsoon_nh = (p_diff.precip > 2.) & (p_clim.precip.sel(month=[5,6,7,8,9]).sum('month')/p_clim.precip.sum('month') > 0.55) 
monsoon_sh = (p_diff.precip < -2.) & (p_clim.precip.sel(month=[11,12,1,2,3]).sum('month')/p_clim.precip.sum('month') > 0.55) 
lats_nh = [p_diff.nlat[i] for i in range(len(p_diff.latitude)) if p_diff.latitude[i] > 0]
lats_sh = [p_diff.nlat[i] for i in range(len(p_diff.latitude)) if p_diff.latitude[i] < 0]

lats_tropics = [p_diff.nlat[i] for i in range(len(p_diff.latitude)) if p_diff.latitude[i] > -30. and p_diff.latitude[i] < 30.]

itcz_aug= np.zeros(len(p_clim.nlon),)
itcz_feb= np.zeros(len(p_clim.nlon),)

for i in range(len(p_clim.nlon)):
    #print(p_clim.latitude[p_clim.precip.sel(month=8)[:,i] == p_clim.precip.sel(month=8)[:,i].max('nlat')])
    itcz_aug[i] = p_clim.latitude[p_clim.precip.sel(month=8)[:,i] == p_clim.precip.sel(month=8, nlat=lats_tropics)[:,i].max('nlat')].values
    itcz_feb[i] = p_clim.latitude[p_clim.precip.sel(month=2)[:,i] == p_clim.precip.sel(month=2, nlat=lats_tropics)[:,i].max('nlat')].values
    
#print(itcz_aug)

# Make plots

# Set figure parameters
rcParams['figure.figsize'] = 12, 7
rcParams['font.size'] = 18

fig = plt.figure(1)
# Contour plot of precipitation difference
f1 = p_diff.precip.plot.contourf(x='longitude', y='latitude', cmap='RdBu', levels=np.arange(-10.,10.1,1.), add_colorbar=False, extend='both')
monsoon_nh.sel(nlat=lats_nh).plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), colors='c', add_colorbar=False, linewidths=2)
monsoon_sh.sel(nlat=lats_sh).plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), colors='c', add_colorbar=False, linewidths=2)
land.lsm[0,:,:].plot.contour(x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k', linewidths=2.)
#plt.plot(p_diff.longitude, itcz_aug, 'r')
#plt.plot(p_diff.longitude, itcz_feb, 'b')
b = plt.quiver(u_diff.lon[::5], u_diff.lat[::2], u_diff.var33.sel(lev=85000.)[::2,::5], v_diff.var34.sel(lev=85000.)[::2,::5], angles='xy', scale=200.)
plt.grid(True,linestyle=':')
plt.ylabel('Latitude')
plt.xlabel('Longitude')
plt.ylim(-60.,60.)
plt.xticks(np.arange(0.,360.,20.))
plt.yticks(np.arange(-60.,61.,30.))

cb1=fig.colorbar(f1, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.15, aspect=30, shrink=0.5)
cb1.set_label('MJJAS - NDJFM precipitation, mm/day')
    
plt.savefig(plot_dir + 'wind_precip_reversal.pdf', format='pdf')
plt.close()

data_u.close()
data_v.close()
data_p.close()