''' 
27/09/2018 Make side by side hm plots of regional precip in GPCP and half_shallow
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from data_handling_updates import month_dic, make_sym

def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons
    
plot_dir = '/scratch/rg419/plots/itcz_plots/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/half_shallow.nc')
data['precipitation'] = make_sym(data.precipitation)*86400.
data_gpcp = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/gpcp_detrendedpentads.nc')
data_gpcp = data_gpcp.rename({'precip_clim': 'precipitation'})

mn_dic = month_dic(1)
tickspace = [13, 31, 49, 68]
labels = ['Mar', 'Jun', 'Sep', 'Dec']    

# Set figure parameters
rcParams['figure.figsize'] = 15, 6
rcParams['font.size'] = 16
    
# Start figure with 6 subplots
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
axes = [ax1, ax2, ax3, ax4, ax5, ax6]
def precip_plot(data, ax, lons, title, levels=np.arange(1.,11.,1.)):
    lons = pick_lons(data,lons)
    f1 = data.precipitation.sel(lon=lons).mean('lon').plot.contourf(ax=ax, x='xofyear', y='lat', levels = levels, add_labels=False, extend='max', cmap='Blues', add_colorbar=False)
    ax.set_ylim(-45,45)
    ax.grid(True,linestyle=':')
    ax.set_yticks(np.arange(-30.,31.,30.))
    ax.set_xlim((1,73))    
    ax.set_xticks(tickspace)
    ax.set_title(title)
    return f1

precip_plot(data_gpcp, ax1, [100.,150.], 'GPCP - East Asia')
precip_plot(data_gpcp, ax2, [70.,100.], 'GPCP - India')
f1 = precip_plot(data_gpcp, ax3, [15.,60.], 'GPCP - South Africa')
precip_plot(data, ax4, [170.,190.], 'Isca - East coast', levels=np.arange(2.,22.,2.))
precip_plot(data, ax5, [80.,100.], 'Isca - Land', levels=np.arange(2.,22.,2.))
f2 = precip_plot(data, ax6, [350.,10.], 'Isca - West coast', levels=np.arange(2.,22.,2.))

ax1.set_ylabel('Latitude')
ax4.set_ylabel('Latitude')
for ax in axes[3:]:
    ax.set_xticklabels(labels,rotation=25)

plt.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.1, hspace=0.25, wspace=0.2)
cb1=fig.colorbar(f1, ax=axes[:3], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.03, aspect=15, shrink=1.)
cb2=fig.colorbar(f2, ax=axes[3:], use_gridspec=True, orientation = 'vertical',fraction=0.05, pad=0.03, aspect=15, shrink=1.)
cb1.set_label('Precipitation, mm/day')
cb2.set_label('Precipitation, mm/day')
    
# Save as a pdf
plt.savefig(plot_dir + 'precip_half_shallow_gpcp.pdf', format='pdf')
plt.close()

data.close()




