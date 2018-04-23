# Plot the surface temperature and precipitation through the year for the monsoon region

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from physics import model_constants as mc, precip_centroid

rcParams['figure.figsize'] = 6, 10
rcParams['font.size'] = 20
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/seasons/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

mn_dic = month_dic(1)

plevels = np.arange(2.,19.,2.)

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=False)
plt.set_cmap('RdBu_r')

tickspace = np.arange(13,72,18)



def precip_panel(run, axisno, period_fac=1.):
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['totp'] = ((data.precipitation)*86400.).mean('lon')
    
    data['mse'] = (data.temp*mc.cp_air + data.sphum*mc.L + data.height*mc.grav).mean('lon')/1000.
    
    p_cent = precip_centroid(run, lat_bound=35.)
    
    f1 = data.totp.plot.contourf(ax=axisno, x='xofyear', y='lat', levels = plevels, add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
    #data.totp.plot.contour(ax=axisno, x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
    p_cent.plot(ax=axisno, color='w')
    cs = data.mse.sel(pfull=850.).plot.contour(ax=axisno, x='xofyear', y='lat', levels=np.arange(200.,401.,20.), add_labels = False, colors='0.7', add_colorbar=False)
    plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
    axisno.set_ylabel('Latitude')
    axisno.set_xlabel('')
    axisno.set_ylim(-35,35)
    axisno.set_yticks(np.arange(-30.,31.,15.))
    axisno.set_xticks(tickspace * period_fac)
    axisno.set_xticklabels('')
    axisno.grid(True,linestyle=':')
    
    
    return f1



f1 = precip_panel('sn_0.500', ax1, period_fac=0.5)    
precip_panel('sn_1.000', ax2)    
precip_panel('sn_2.000', ax3, period_fac = 2.)    
    
labels = [mn_dic[(int(k)+5)/6 ] for k in tickspace]
ax3.set_xticks(tickspace * 2.)
ax3.set_xticklabels(labels,rotation=25)


plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=fig.colorbar(f1, ax=(ax1, ax2, ax3), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
cb1.set_label('Precipitation, mm/day')
 
plt.savefig(plot_dir+'precip_mse_hm_centroid.pdf', format='pdf')
plt.close()


# Also plot centroid in 'real time' to allow comparison

# Identify when centroid is over the Equator:
def get_pcent_eq(run):
    p_cent = precip_centroid(run, lat_bound=35.)

    eq_time = int(np.abs(p_cent[20:]).argmin('xofyear') + 20)
    pcent_eq = p_cent[eq_time-5:eq_time+15]
    return pcent_eq
    
p_cent_05 = precip_centroid('sn_0.500', lat_bound=35.)
p_cent_10 = precip_centroid('sn_1.000', lat_bound=35.)
p_cent_20 = precip_centroid('sn_2.000', lat_bound=35.)

rcParams['figure.figsize'] = 9, 4

plt.plot(p_cent_05[20:36], 'b', LineWidth=2)
plt.plot(p_cent_10[39:55], 'k', LineWidth=2)
plt.plot(p_cent_20[78:94], 'r', LineWidth=2)
#plt.plot([20,20], [-30,30], '--k')
plt.xlim([0,15])
plt.ylabel('Latitude')
plt.xlabel('Time')
plt.tight_layout()
plt.savefig(plot_dir+'real_time_centroid.pdf', format='pdf')
plt.close()

