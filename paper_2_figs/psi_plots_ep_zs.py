"""
Plot precip hms for eddy permitting and zonally symmetric runs (26/03/2018)

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_mse_plot
from pylab import rcParams
from hadley_cell import mass_streamfunction
from data_handling_updates import model_constants as mc

data_sn1 = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/sn_1.000.nc')
data_sn1_zs = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/sn_1_sst_zs.nc')
data_sine = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/sine_sst_10m.nc')
data_sine_zs = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/sine_sst_10m_zs.nc')

    
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
#rcParams['figure.figsize'] = 10, 10
rcParams['figure.figsize'] = 10, 6
rcParams['font.size'] = 14


def psi_u_plot(data, tf, ax):
    
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    
    n = len(data.xofyear.values)//2
    psi_temp = np.zeros(psi.values.shape)
    for i in range(0,n):
        psi_temp[:,i,:] = (psi[:,i,:].values - psi[::-1,i+n,:].values)/2.
        psi_temp[:,i+n,:] = -1.*psi_temp[::-1,i,:]
    psi = xr.DataArray(psi_temp, coords=[data.lat, data.xofyear.values, psi.pfull], dims=['lat', 'xofyear', 'pfull'])

    m = mc.omega * mc.a**2. * np.cos(psi.lat*np.pi/180.)**2. + data.ucomp.mean('lon') * mc.a * np.cos(psi.lat*np.pi/180.)
    
    m_levs = mc.omega * mc.a**2. * np.cos(np.arange(-60.,1.,5.)*np.pi/180.)**2.
    
    f1 = data.ucomp[tf[0]:tf[1],:,:].mean(('xofyear','lon')).plot.contourf(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-50.,50.1,5.), extend='both', add_labels=False, add_colorbar=False)
    
    m[:,tf[0]:tf[1],:].mean('xofyear').plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=m_levs, colors='0.7', add_labels=False)
    
    psi[:,tf[0]:tf[1],0:19].mean('xofyear').plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(0.,601,100.), colors='k', add_labels=False)
    psi[:,tf[0]:tf[1],0:19].mean('xofyear').plot.contour(ax=ax, x='lat', y='pfull', yincrease=False, levels=np.arange(-600.,0.,100.), colors='k', linestyles='dashed', add_labels=False)
        
    ax.set_xlim(-35,35)
    ax.set_xticks(np.arange(-30,31,15))
    ax.grid(True,linestyle=':')
    
    return f1


#fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9), (ax10, ax11, ax12)) = plt.subplots(4, 3, sharex='col', sharey='row')
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
plt.set_cmap('RdBu_r')
    

f1=psi_u_plot(data_sn1, [31,35], ax=ax1)
psi_u_plot(data_sn1, [38,42], ax=ax2)
psi_u_plot(data_sn1, [45,49], ax=ax3)

#psi_u_plot(data_sine, [31,35], ax=ax4)
#psi_u_plot(data_sine, [38,42], ax=ax5)
#psi_u_plot(data_sine, [45,49], ax=ax6)

psi_u_plot(data_sn1_zs, [31,35], ax=ax4)
psi_u_plot(data_sn1_zs, [38,42], ax=ax5)
psi_u_plot(data_sn1_zs, [45,49], ax=ax6)

#psi_u_plot(data_sine_zs, [31,35], ax=ax10)
#psi_u_plot(data_sine_zs, [38,42], ax=ax11)
#psi_u_plot(data_sine_zs, [45,49], ax=ax12)

#for ax in [ax10, ax11, ax12]:
for ax in [ax4, ax5, ax6]:
    ax.set_xlabel('Latitude')
    
for ax in [ax1, ax4]: #, ax7, ax10]:
    ax.set_ylabel('Pressure, hPa')
    ax.set_yticks([0,200,400,600,800,1000.])
    
ax1.set_title('Pentad 32-35', fontsize=14)
ax2.set_title('Pentad 39-42', fontsize=14)
ax3.set_title('Pentad 46-49', fontsize=14)

i=0
labels=['b)','c)','e)','f)']
for ax in [ax2, ax3, ax5, ax6]:
    ax.text(-45, 0., labels[i])
    i=i+1
ax1.text(-50, 0., 'a)')
ax4.text(-50, 0., 'd)')


#axes = (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12)
axes = (ax1, ax2, ax3, ax4, ax5, ax6)


plt.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.03, hspace=0.2, wspace=0.2)
#Colorbar
cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.1, pad=0.12, aspect=30, shrink=0.5)
cb1.set_label('Zonal wind speed, m/s')

plt.savefig(plot_dir+'psi_ep_zs.pdf', format='pdf')
plt.close()        
