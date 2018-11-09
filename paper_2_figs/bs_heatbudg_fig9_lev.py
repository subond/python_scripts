"""
Reproduce figures 9 from SB08 using data from the updated sn_1.000 run at a given level (01/09/2018)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, model_constants as mc
from climatology import precip_centroid
from pylab import rcParams
from hadley_cell import mass_streamfunction

lev=850.
plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
rcParams['figure.figsize'] = 5.5, 7.
rcParams['font.size'] = 14

def fig_9(run, ax, pentad=40, lev=850.):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    convTtotheta=(1000./data.pfull)**mc.kappa
    
    heating = data.dt_tg_condensation + data.dt_tg_convection + data.dt_tg_diffusion + data.tdt_rad
    heating_theta = (heating*convTtotheta).mean('lon')*86400.
    
    theta = data.temp*convTtotheta
    dthetady = gr.ddy(theta.mean('lon'), vector=False)
    dthetadp = gr.ddp(theta.mean('lon'))
    dthetadt = gr.ddt(theta.mean('lon'))*86400.
    
    vcomp_theta = data.vcomp_temp*convTtotheta
    vcomp_theta_eddy = vcomp_theta.mean('lon') - data.vcomp.mean('lon')*theta.mean('lon')
    
    omega_theta = data.omega_temp*convTtotheta
    omega_theta_eddy = omega_theta.mean('lon') - data.omega.mean('lon')*theta.mean('lon')
    
    vdthetady_mean = -1.*data.vcomp.mean('lon') * dthetady*86400.
    wdthetadp_mean = -1.*data.omega.mean('lon') * dthetadp*86400.
    
    div_vt_eddy = -1.*gr.ddy(vcomp_theta_eddy, vector=True)*86400.
    div_wt_eddy = -1.*gr.ddp(omega_theta_eddy)*86400.
    
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -40. and data.lat[i] <= 40.]
    
    #heating_theta.sel(pfull=lev, xofyear=pentad, lat=lats).plot(ax=ax, color='k')
    dthetadt.sel(pfull=lev, xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    vdthetady_mean.sel(pfull=lev, xofyear=pentad, lat=lats).plot(ax=ax, color='C1')
    (wdthetadp_mean + heating_theta).sel(pfull=lev, xofyear=pentad, lat=lats).plot(ax=ax, color='C2')
   # (vdthetady_mean + wdthetadp_mean).sel(pfull=lev, xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    #(-dthetadt + vdthetady_mean + wdthetadp_mean + heating_theta + div_vt_eddy + div_wt_eddy).sel(pfull=lev, xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    #(div_vt_eddy + div_wt_eddy).sel(pfull=lev, xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylim(-0.25,0.25)
    #ax.set_ylim(-5.,5.)
    ax.set_xlim(-30.,30.)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('Heating rate, K/day')
    

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)

fig_9('sn_1.000_evap_fluxes_heattrans', ax=ax1, pentad=35, lev=lev)
fig_9('sn_1.000_evap_fluxes_heattrans', ax=ax2, pentad=40, lev=lev)
fig_9('sn_1.000_evap_fluxes_heattrans', ax=ax3, pentad=45, lev=lev)
fig_9('sn_1.000_evap_fluxes_heattrans', ax=ax4, pentad=50, lev=lev)

ax4.set_xlabel('Latitude')

ax1.text(-55, 500., 'a)')
ax2.text(-55, 500., 'b)')
ax3.text(-55, 500., 'c)')
ax4.text(-55, 500., 'd)')

plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)

plt.savefig(plot_dir+'sb08_fig9_lev_' + str(int(lev)) + '.pdf', format='pdf')
plt.close()
