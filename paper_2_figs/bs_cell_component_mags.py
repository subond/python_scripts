"""
Reproduce figure 9 from SB10 using data from control and zs runs (2/1/2019)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, model_constants as mc, make_sym
from climatology import precip_centroid
from hadley_cell import get_edge_psi
from pylab import rcParams
from hadley_cell import mass_streamfunction
from pcent_rate_max import p_cent_rate_max

        
def fig_9(run, ax, rotfac=1., lev=200., check=False, eddies=True):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    psi_full = mass_streamfunction(data, a=6376.0e3, dp_in=50.) / 1.e9
    vcomp = data.vcomp
    
    sinphi = np.sin(data.lat * np.pi/180.)
    cosphi = np.cos(data.lat * np.pi/180.)
    coriolis = 2.* rotfac * mc.omega * sinphi 
    
    dudt = gr.ddt(data.ucomp)
    dudy = gr.ddy(data.ucomp)
    duvdy = gr.ddy(data.ucomp_vcomp, uv=True)
    dudp = gr.ddp(data.ucomp)
    duwdp = gr.ddp(data.ucomp_omega)
    
    vdudy_mean = data.vcomp * dudy
    wdudp_mean = data.omega * dudp
    adv_v = (vdudy_mean + wdudp_mean)/coriolis
    eddy_v = (duvdy + duwdp)/coriolis - adv_v
    time_v = dudt/coriolis
        
    def psi_comp(comp, v):
        data['vcomp'] = comp
        psi_comp = mass_streamfunction(data, a=6376.0e3, dp_in=50.) / 1.e9
        data['vcomp'] = v
        return psi_comp
    
    psi_eddy = psi_comp(eddy_v, vcomp)
    psi_adv = psi_comp(adv_v, vcomp)
    psi_time = psi_comp(time_v, vcomp)
    
    psi_max = np.zeros(72)
    psi_max_loc = np.zeros(72)
    psi_eddy_max = np.zeros(72)
    psi_adv_max = np.zeros(72)
    psi_time_max = np.zeros(72)
    
    pad = 3
    i1 = 32 - pad #29
    i2 = 32 + pad
    
    for i in range(1,73):
        psi_full_i = psi_full.sel(xofyear=i, pfull=lev)
        psi_full_i[i1:i2] = 0. # set to 0 for lats near equator: eqward of -7:7 degrees
        psi_eddy[i1:i2] = 0.
        psi_adv[i1:i2] = 0.
        psi_time[i1:i2] = 0.
        psi_max[i-1] = -1.* psi_full_i.where(psi_full_i == psi_full_i.min(), drop=True).values
        psi_max_loc[i-1] = psi_full_i.lat.where(psi_full_i == psi_full_i.min(), drop=True)
        psi_eddy_max[i-1] = -1.* psi_eddy.sel(xofyear=i, pfull=lev).where(psi_full_i == psi_full_i.min(), drop=True).values
        psi_adv_max[i-1]  = -1.* psi_adv.sel(xofyear=i, pfull=lev).where(psi_full_i == psi_full_i.min(), drop=True).values
        psi_time_max[i-1] = -1.* psi_time.sel(xofyear=i, pfull=lev).where(psi_full_i == psi_full_i.min(), drop=True).values
    
    if check:
        psi_full.sel(pfull=lev).plot.contourf(x='xofyear', y='lat', levels=np.arange(-400.,401.,50.))
        plt.plot(np.arange(1,73), psi_max_loc, color='k')
        plt.xticks(np.arange(0., 73., 12.))
        
        plt.figure(2)
        psi_adv.sel(pfull=lev).plot.contourf(x='xofyear', y='lat', levels=np.arange(-400.,401.,50.))
        plt.plot(np.arange(1,73), psi_max_loc, color='k')
        plt.xticks(np.arange(0., 73.,12.))
        
        plt.figure(3)
        psi_eddy.sel(pfull=lev).plot.contourf(x='xofyear', y='lat', levels=np.arange(-400.,401.,50.))
        plt.plot(np.arange(1,73), psi_max_loc, color='k')
        plt.xticks(np.arange(0., 73.,12.))
        
        plt.show()
    
    ax.plot(np.arange(1,73), psi_max, color='k', ls ='--')
    ax.plot(np.arange(1,73), psi_adv_max, color='k')
    ax.plot(np.arange(1,73), psi_time_max, color='0.7')
    if eddies:
        ax.plot(np.arange(1,73), psi_eddy_max, color='k', ls=':')
    ax.set_xticks(np.arange(0., 73.,12.))
    ax.set_xlim(0.,72.)


# Set plotting directory
plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

# Set figure parameters
rcParams['figure.figsize'] = 5.5, 7.5
rcParams['font.size'] = 14

fig, ((ax1), (ax2)) = plt.subplots(2, 1)


fig_9('sn_1.000_evap_fluxes_heattrans', ax=ax1)
fig_9('sn_1_sst_zs', ax=ax2, eddies=False)

ax1.text(-5, 350., 'a)')
ax2.text(-5, 660., 'b)')
ax1.legend(['$\Psi$', '$\Psi_M$', '$\Psi_t$', '$\Psi_E$'])
ax2.set_xlabel('Pentad')
ax1.set_ylabel('$\Psi$, (10$^9$kgs$^{-1}$)')
ax2.set_ylabel('$\Psi$, (10$^9$kgs$^{-1}$)')
ax1.set_title('control')
ax2.set_title('control-zs')

plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.07, hspace=0.3)

plt.savefig(plot_dir+'bs10_fig9.pdf', format='pdf')
plt.close()
plt.show()
