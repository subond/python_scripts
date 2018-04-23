"""
Create figure showing pressure mean v at upper and lower levels at 150 hPa and 850 hPa for (top) 10 N, (bottom) 40 S during monsoon period.

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from pylab import rcParams

    
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
        
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/full_qflux.nc')

lat_hm = data.lat[np.argmin(np.abs(data.lat - 10.))]

vmass_10_lower_b = data.vcomp.sel(lat=lat_hm)[18:18+4,0:9,:].mean(('xofyear','pfull'))
vmass_10_upper_b = -1.*data.vcomp.sel(lat=lat_hm)[18:18+4,10:19,:].mean(('xofyear','pfull'))

# Get correlation coefficient
lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= 60. and data.lon[i] < 150.]
print np.corrcoef(vmass_10_lower_b.sel(lon=lons), vmass_10_upper_b.sel(lon=lons))[0,1]**2
print np.corrcoef(vmass_10_lower_b, vmass_10_upper_b)[0,1]**2

vmass_10_lower_a = data.vcomp.sel(lat=lat_hm)[39:39+4,0:9,:].mean(('xofyear','pfull'))
vmass_10_upper_a = -1.*data.vcomp.sel(lat=lat_hm)[39:39+4,10:19,:].mean(('xofyear','pfull'))

print np.corrcoef(vmass_10_lower_a.sel(lon=lons), vmass_10_upper_a.sel(lon=lons))[0,1]**2
print np.corrcoef(vmass_10_lower_a, vmass_10_upper_a)[0,1]**2


lat_hm = data.lat[np.argmin(np.abs(data.lat - 30.))]
vmass_40_lower = data.vcomp.sel(lat=lat_hm)[39:39+4,0:9,:].mean(('xofyear','pfull'))
vmass_40_upper = -1.*data.vcomp.sel(lat=lat_hm)[39:39+4,10:19,:].mean(('xofyear','pfull'))

print np.corrcoef(vmass_40_lower.sel(lon=lons), vmass_40_upper.sel(lon=lons))[0,1]**2
print np.corrcoef(vmass_40_lower, vmass_40_upper)[0,1]**2


rcParams['figure.figsize'] = 6, 7
# Two subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
#First plot
ax1.fill_between([60.,150.], -6., 2., facecolor='gray', alpha=0.3, edgecolor="none")
vmass_10_lower_b.plot(ax=ax1, color='r', linewidth=2)
vmass_10_upper_b.plot(ax=ax1, color='b', linewidth=2)
ax1.set_ylabel('$v$, m/s')
#ax1.set_ylim([-6,6])
ax1.set_xlim([0,360])
ax1.set_title('')
ax1.set_xlabel('')
ax1.set_yticks(np.arange(-6,2.1,2.))
ax1.grid(True,linestyle=':')
ax1.text(-85, 2, 'a)')

ax2.fill_between([60.,150.], -2., 6., facecolor='gray', alpha=0.3, edgecolor="none")
vmass_10_lower_a.plot(ax=ax2, color='r', linewidth=2)
vmass_10_upper_a.plot(ax=ax2, color='b', linewidth=2)
ax2.set_ylabel('$v$, m/s')
#ax2.set_ylim([-6,6])
ax2.set_xlim([0,360])
ax2.set_title('')
ax2.set_xlabel('')
ax2.set_yticks(np.arange(-2,6.1,2.))
ax2.grid(True,linestyle=':')
ax2.text(-85, 6, 'b)')

ax3.fill_between([60.,150.], -8., 8., facecolor='gray', alpha=0.3, edgecolor="none")
vmass_40_lower.plot(ax=ax3, color='r', linewidth=2)
vmass_40_upper.plot(ax=ax3, color='b', linewidth=2)
ax3.set_ylabel('$v$, m/s')
#ax3.set_ylim([-6,6])
ax3.set_xlim([0,360])
ax3.set_title('')
ax3.set_xlabel('Longitude')
ax3.set_yticks(np.arange(-8,8.1,4.))
ax3.set_xticks(np.arange(0,360.1,60.))
ax3.grid(True,linestyle=':')
ax3.text(-85, 8, 'c)')

plt.subplots_adjust(left=0.17, right=0.83, top=0.95, bottom=0.1, hspace=0.2)
plt.savefig(plot_dir+'vmass_figure.pdf', format='pdf')
plt.close()
    


