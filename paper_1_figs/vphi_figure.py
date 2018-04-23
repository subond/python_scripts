"""
Create figure showing mean geopotential advection, and contribution from mean state to this, at 150 hPa and 850 hPa for (left) 10 N, (right) 40 S during monsoon period.

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from pylab import rcParams

    
rcParams['figure.figsize'] = 6, 9
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
        
data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/geopot_flux_full_qflux.nc')
vphi_mean = data.phi * data.vcomp

lat_hm = data.lat[np.argmin(np.abs(data.lat - 10.))]

vphi_150_10 = (data.vphi.sel(pfull=150.,lat=lat_hm)[39:39+4,:].mean('pentad'))/10.**4
vphi_850_10 = (data.vphi.sel(pfull=850.,lat=lat_hm)[39:39+4,:].mean('pentad'))/10.**4

vphi_mean_150_10 = (vphi_mean.sel(pfull=150.,lat=lat_hm)[39:39+4,:].mean('pentad'))/10.**4
vphi_mean_850_10 = (vphi_mean.sel(pfull=850.,lat=lat_hm)[39:39+4,:].mean('pentad'))/10.**4

lat_hm = data.lat[np.argmin(np.abs(data.lat + 40.))]

vphi_150_40 = (data.vphi.sel(pfull=150.,lat=lat_hm)[39:39+4,:].mean('pentad'))/10.**4
vphi_850_40 = (data.vphi.sel(pfull=850.,lat=lat_hm)[39:39+4,:].mean('pentad'))/10.**4

vphi_mean_150_40 = (vphi_mean.sel(pfull=150.,lat=lat_hm)[39:39+4,:].mean('pentad'))/10.**4
vphi_mean_850_40 = (vphi_mean.sel(pfull=850.,lat=lat_hm)[39:39+4,:].mean('pentad'))/10.**4



# Two subplots
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
#First plot
vphi_150_10.plot(ax=ax1, color='r', linewidth=2)
#vphi_mean_150_10.plot(ax=ax1, color='k', ls='--', linewidth=2)
ax1.set_ylabel('10$^4$ v$\Phi$, 150 hPa', color='r')
ax1.set_ylim([-240,240])
ax1.set_xlim([0,360])
ax1.set_title('')
ax1.set_xlabel('')
ax1.set_yticks(np.arange(-240.,241.,80.))
for tl in ax1.get_yticklabels():
    tl.set_color('r')

ax3 = ax1.twinx()
vphi_850_10.plot(ax=ax3, color='b', linewidth=2)
#vphi_mean_850_10.plot(ax=ax3, color='k', ls='--', linewidth=2)
ax3.set_ylabel('10$^4$ v$\Phi$, 850 hPa', color='b')
ax3.set_ylim([-15,15])
ax3.set_xlim([0,360])
ax3.set_title('')
for tl in ax3.get_yticklabels():
    tl.set_color('b')
    
ax1.grid(True,linestyle=':')
ax1.text(-50, 250, 'a)')


vphi_150_40.plot(ax=ax2, color='r', linewidth=2)
#vphi_mean_150_40.plot(ax=ax2, color='k', ls='--', linewidth=2)
ax2.set_ylabel('10$^4$ v$\Phi$, 150 hPa', color='r')
ax2.set_ylim([-50,50])
ax2.set_xlim([0,360])
ax2.set_title('')
ax2.set_xlabel('Longitude')
#ax2.set_yticks(np.arange(-240.,241.,80.))
for tl in ax2.get_yticklabels():
    tl.set_color('r')

ax4 = ax2.twinx()
vphi_850_40.plot(ax=ax4, color='b', linewidth=2)
#vphi_mean_850_40.plot(ax=ax4, color='k', ls='--', linewidth=2)
ax4.set_ylabel('10$^4$ v$\Phi$, 850 hPa', color='b')
ax4.set_ylim([-5,5])
ax4.set_xlim([0,360])
ax4.set_title('')
for tl in ax4.get_yticklabels():
    tl.set_color('b')
    
ax2.grid(True,linestyle=':')
ax2.text(-50, 250, 'b)')

plt.subplots_adjust(left=0.17, right=0.83, top=0.95, bottom=0.1, hspace=0.1)


plt.savefig(plot_dir+'vphi_figure.pdf', format='pdf')
plt.close()
    


