''' 
11/10/2018 Make CMAP hms of different monsoon regions:
- East Asia
- Western Pacific
- India
- West Africa
- South Africa
- South America
- North America
- Bay of Bengal
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

data = xr.open_dataset('/disca/share/rg419/CMAP_precip.pentad.mean.nc')
data.coords['pentad'] = (('time'), np.tile(np.arange(1,74),38))
data = data.groupby('pentad').mean('time')
    
mn_dic = month_dic(1)
tickspace_nh = [1, 19, 37, 55]
labels_nh = ['1st Jan', '1st Apr', '30th Jun', '28th Sept']    

tickspace_sh = [-36, -18, 1, 19]
labels_sh = ['30th Jun', '28th Sept', '1st Jan', '1st Apr']    

# Set figure parameters
rcParams['figure.figsize'] = 15, 6
rcParams['font.size'] = 14
    
# Start figure with 6 subplots
fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex='col', sharey='row')
axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
def precip_plot(data, ax, lons, title, levels=np.arange(0.,15.,2.), sh=False):
    lons = pick_lons(data,lons)
    if sh:
        precip = xr.DataArray(data.precip.sel(lon=lons).mean('lon').values, [('pentad', (np.arange(1,74)-37)%73-36), ('lat', data.lat)])
        precip = precip.sortby('pentad')
        f1 = precip.plot.contourf(ax=ax, x='pentad', y='lat', levels = levels, add_labels=False, extend='max', cmap='Blues', add_colorbar=False)
    else:
         f1 = data.precip.sel(lon=lons).mean('lon').plot.contourf(ax=ax, x='pentad', y='lat', levels = levels, add_labels=False, extend='max', cmap='Blues', add_colorbar=False)
    ax.set_ylim(-45,45)
    ax.grid(True,linestyle=':')
    ax.set_yticks(np.arange(-30.,31.,30.))
    ax.set_title(title)
    return f1

f1 = precip_plot(data, ax1, [60.,80.], 'India')
precip_plot(data, ax2, [90.,100.], 'Bay of Bengal')
precip_plot(data, ax3, [110.,140.], 'East Asia')
precip_plot(data, ax4, [15.,50.], 'South Africa', sh=True)
precip_plot(data, ax5, [140.,160.], 'Western North Pacific')
precip_plot(data, ax6, [240.,280.], 'North America')
precip_plot(data, ax7, [340.,20.], 'West Africa')
precip_plot(data, ax8, [280.,320.], 'South America', sh=True)


ax1.set_ylabel('Latitude')
ax5.set_ylabel('Latitude')
for ax in [ax1,ax2,ax3,ax5,ax6,ax7]:
    ax.set_xticks(tickspace_nh)
    ax.set_xticklabels(labels_nh,rotation=25)
for ax in [ax4,ax8]:
    ax.set_xticks(tickspace_sh)
    ax.set_xticklabels(labels_sh,rotation=25)

plt.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.2)
cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.15, aspect=60, shrink=0.5)
cb1.set_label('Precipitation, mm/day')
    
# Save as a pdf
plt.savefig(plot_dir + 'precip_monsoons_cmap.pdf', format='pdf')
plt.close()

data.close()




