# Plot total vorticity tendency (lat-time) - use time gradient of vorticity to do so.

from data_handling import month_dic
from physics import gradients as gr
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh

rcParams['figure.figsize'] = 6, 10
rcParams['font.size'] = 20
rcParams['text.usetex'] = True
    

plot_dir = '/scratch/rg419/plots/seasons/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

lev=150.

data_360 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_sn_1.000.nc')
data_720 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_sn_2.000.nc')
data_180 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_sn_0.500.nc')

def get_abs_vort(data):
    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    data['abs_vort'] = vor.mean('lon')*86400.
    data['vort_tend'] = gr.ddt(data.abs_vort, timedir='pentad')*86400.

get_abs_vort(data_180)
get_abs_vort(data_360)
get_abs_vort(data_720)

mn_dic = month_dic(1)
tickspace = np.arange(13,72,18)
labels = [mn_dic[(int(k)+5)/6 ] for k in tickspace]    

levels = np.arange(-0.1,0.11,0.02)


# Begin plotting: 3 subplots

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=False)
plt.set_cmap('RdBu_r')
#First plot
f1 = data_180.vort_tend.sel(pfull=lev).plot.contourf(ax=ax1, x='pentad', y='lat', levels = levels, add_colorbar=False, add_labels=False, extend='both')
ax1.contour(data_180.pentad, data_180.lat, data_180.abs_vort.sel(pfull=lev).T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
ax1.set_ylabel('Latitude')
ax1.set_ylim(-60,60)
ax1.set_yticks(np.arange(-60.,61.,30.))
ax1.grid(True,linestyle=':')
ax1.set_xticks(tickspace*0.5)
ax1.set_xticklabels('')

f1 = data_360.vort_tend.sel(pfull=lev).plot.contourf(ax=ax2, x='pentad', y='lat', levels = levels, add_colorbar=False, add_labels=False, extend='both')
ax2.contour(data_360.pentad, data_360.lat, data_360.abs_vort.sel(pfull=lev).T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
ax2.set_ylabel('Latitude')
ax2.set_ylim(-60,60)
ax2.set_yticks(np.arange(-60.,61.,30.))
ax2.grid(True,linestyle=':')
ax2.set_xticks(tickspace)
ax2.set_xticklabels('')

f1 = data_720.vort_tend.sel(pfull=lev).plot.contourf(ax=ax3, x='pentad', y='lat', levels = levels, add_colorbar=False, add_labels=False, extend='both')
ax3.contour(data_720.pentad, data_720.lat, data_720.abs_vort.sel(pfull=lev).T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
ax3.set_ylabel('Latitude')
ax3.set_ylim(-60,60)
ax3.set_yticks(np.arange(-60.,61.,30.))
ax3.grid(True,linestyle=':')

#ax3.set_xlim((1,72))
ax3.set_xlabel('')
ax3.set_xticks(tickspace*2.)
ax3.set_xticklabels(labels,rotation=25)
    
plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
#Colorbar
cb1=fig.colorbar(f1, ax=(ax1, ax2, ax3), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
cb1.set_label('Vorticity tendency, day$^{-2}$')
 
plt.savefig(plot_dir+'vort_total_alt.pdf', format='pdf')
plt.close()        

