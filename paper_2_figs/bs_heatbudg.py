"""
Reproduce figures 8 and 9 from SB08 using data from the updated sn_1.000 run (01/09/2018)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, model_constants as mc
from climatology import precip_centroid
from pylab import rcParams
from hadley_cell import mass_streamfunction


def fig_8(run, ax, pentad=45):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    heating = data.dt_tg_condensation + data.dt_tg_convection + data.dt_tg_diffusion + data.tdt_rad
    
    convTtotheta=(1000./data.pfull)**mc.kappa
    theta = data.temp*convTtotheta
    heating_theta = heating*convTtotheta
    
    dthetady = gr.ddy(theta.mean('lon'), vector=False)
    dthetadp = gr.ddp(theta.mean('lon'))
    dthetadt = gr.ddt(theta.mean('lon'))
    vdthetady_mean = data.vcomp.mean('lon') * dthetady
    wdthetadp_mean = data.omega.mean('lon') * dthetadp
    adv_heating = -1. *(vdthetady_mean + wdthetadp_mean)
    
    sinphi = np.sin(data.lat * np.pi/180.)
    cosphi = np.cos(data.lat * np.pi/180.)
    
    theta_850 = theta.sel(pfull=850.).mean('lon')
    theta_max = theta_850.max('lat')
    theta_max_lat = data.lat.values[theta_850.argmax('lat').values]
    theta_max_lat = xr.DataArray(theta_max_lat, coords=[data.xofyear], dims=['xofyear'])
    sinphimax = np.sin(theta_max_lat * np.pi/180.)
    
    theta_m = theta_max - 300.*mc.omega**2.*mc.a**2./(2.*mc.grav*14000.) * (sinphi**2.-sinphimax**2.)**2./cosphi**2.
    
    theta_rc = theta_850 + heating_theta.sel(pfull=850.).mean('lon') * 86400.

    theta_adv = theta_850 + adv_heating.sel(pfull=850.) * 86400. 

    theta_net = theta_850 + dthetadt.sel(pfull=850.) * 86400. * 10.
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -30. and data.lat[i] <= 30.]
    
    theta_850.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k')
    theta_m.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    theta_rc.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    theta_adv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='r', linestyle='-.')
    theta_net.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C1', linestyle='--')
    ax.set_title('')
    ax.set_xlabel('')
    #ax.set_ylim(285.,315.)
    ax.set_ylim(290.,310.)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('$\Theta$, K')


plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
rcParams['figure.figsize'] = 5.5, 7.
rcParams['font.size'] = 14
    
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)

fig_8('sn_1.000_evap_fluxes_heattrans', ax=ax1, pentad=35)
fig_8('sn_1.000_evap_fluxes_heattrans', ax=ax2, pentad=40)
fig_8('sn_1.000_evap_fluxes_heattrans', ax=ax3, pentad=45)
fig_8('sn_1.000_evap_fluxes_heattrans', ax=ax4, pentad=50)

ax4.set_xlabel('Latitude')

ax1.text(-55, 315., 'a)')
ax2.text(-55, 315., 'b)')
ax3.text(-55, 315., 'c)')
ax4.text(-55, 315., 'd)')

plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)

plt.savefig(plot_dir+'sb08_fig8.pdf', format='pdf')
plt.close()
    


def fig_9(run, ax, pentad=40):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    convTtotheta=(1000./data.pfull)**mc.kappa
    
    heating = data.dt_tg_condensation + data.dt_tg_convection + data.dt_tg_diffusion + data.tdt_rad
    rho = data.pfull*100./mc.rdgas/data.temp
    heating_theta = (heating*convTtotheta).mean('lon')
    
    theta = data.temp*convTtotheta
    dthetady = gr.ddy(theta.mean('lon'), vector=False)
    dthetadp = gr.ddp(theta.mean('lon'))
    dthetadt = gr.ddt(theta.mean('lon'))
    
    vcomp_theta = data.vcomp_temp*convTtotheta
    vcomp_theta_eddy = vcomp_theta.mean('lon') - data.vcomp.mean('lon')*theta.mean('lon')
    
    vdthetady_mean = data.vcomp.mean('lon') * dthetady
    wdthetadp_mean = data.omega.mean('lon') * dthetadp
    
    def column_int(var_in):
        var_int = mc.cp_air * var_in.sum('pfull')*5000./mc.grav
        return var_int
    
    vdtdy_mean_int = -1. * column_int(vdthetady_mean)
    wdtdp_mean_int = -1. * column_int(wdthetadp_mean)
    
    vt_eddy_int = -1. * column_int(vcomp_theta_eddy)
    div_vt_eddy_int = gr.ddy(vt_eddy_int, vector=True)
    
    
    heating_theta_int = column_int(heating_theta)
    dthetadt_int = column_int(dthetadt)
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -40. and data.lat[i] <= 40.]
    
    heating_theta_int.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k')
    dthetadt_int.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='r')
    #wdtdp_mean_int.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C1')
    #vdtdy_mean_int.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    (vdtdy_mean_int + wdtdp_mean_int).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    #(-dthetadt_int + vdtdy_mean_int + wdtdp_mean_int + heating_theta_int + div_vt_eddy_int).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    (vdtdy_mean_int + wdtdp_mean_int + heating_theta_int + div_vt_eddy_int).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    div_vt_eddy_int.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylim(-500.,500.)
    ax.set_xlim(-30.,30.)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('Heating rate, W/m$^2$')
    

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)

fig_9('sn_1.000_evap_fluxes_heattrans', ax=ax1, pentad=35)
fig_9('sn_1.000_evap_fluxes_heattrans', ax=ax2, pentad=40)
fig_9('sn_1.000_evap_fluxes_heattrans', ax=ax3, pentad=45)
fig_9('sn_1.000_evap_fluxes_heattrans', ax=ax4, pentad=50)

ax4.set_xlabel('Latitude')

ax1.text(-55, 500., 'a)')
ax2.text(-55, 500., 'b)')
ax3.text(-55, 500., 'c)')
ax4.text(-55, 500., 'd)')

plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)

plt.savefig(plot_dir+'sb08_fig9.pdf', format='pdf')
plt.close()
