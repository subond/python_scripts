"""
Reproduce moist versions of figures 8 and 9 from SB08 using data from the updated ap_2 run (06/09/2018)
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
    hum_tend = data.dt_qg_condensation + data.dt_qg_convection + data.dt_qg_diffusion
    
    convTtotheta=(1000./data.pfull)**mc.kappa
    theta = data.temp*convTtotheta
    heating_theta = heating*convTtotheta
    heating_equiv_theta = heating_theta + mc.L/mc.cp_air * hum_tend * convTtotheta
    
    sinphi = np.sin(data.lat * np.pi/180.)
    cosphi = np.cos(data.lat * np.pi/180.)
    
    theta_850 = theta.sel(pfull=850.).mean('lon')
    
    theta_equiv = (data.temp + mc.L/mc.cp_air * data.sphum/(1-data.sphum)) * convTtotheta
    theta_equiv_850 = theta_equiv.sel(pfull=850.).mean('lon')
    
    dthetady_equiv = gr.ddy(theta_equiv.mean('lon'), vector=False)
    dthetadp_equiv = gr.ddp(theta_equiv.mean('lon'))
    dthetadt_equiv = gr.ddt(theta_equiv.mean('lon'))
    vdthetady_mean = data.vcomp.mean('lon') * dthetady_equiv
    wdthetadp_mean = data.omega.mean('lon') * dthetadp_equiv
    adv_heating_equiv = -1. *(vdthetady_mean + wdthetadp_mean)
    
    def get_latmax(data_in):
        data_max = data_in.max('lat')
        data_max_lat = data_in.lat.values[data_in.argmax('lat').values]
        data_max_lat = xr.DataArray(data_max_lat, coords=[data_in.xofyear], dims=['xofyear'])
        return data_max, data_max_lat
    
    theta_max, theta_max_lat = get_latmax(theta_850)
    #theta_equiv_max, theta_equiv_max_lat = get_latmax(theta_equiv_850)
    sinphimax = np.sin(theta_max_lat * np.pi/180.)
    #sinphimax_equiv = np.sin(theta_equiv_max_lat * np.pi/180.)
    theta_m = theta_max - 300.*mc.omega**2.*mc.a**2./(2.*mc.grav*14000.) * (sinphi**2.-sinphimax**2.)**2./cosphi**2.
    #theta_m_equiv = theta_equiv_max - 300.*mc.omega**2.*mc.a**2./(2.*mc.grav*14000.) * (sinphi**2.-sinphimax_equiv**2.)**2./cosphi**2.
    
    theta_rc = theta_850 + heating_theta.sel(pfull=850.).mean('lon') * 86400. *5.
    
    theta_equiv_rc = theta_equiv_850 + heating_equiv_theta.sel(pfull=850).mean('lon') * 86400. *5.
    
    theta_equiv_adv = theta_equiv_850 + adv_heating_equiv.sel(pfull=850.) * 86400. *5.

    theta_equiv_net = theta_equiv_850 + dthetadt_equiv.sel(pfull=850.) * 86400.*5.
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -60. and data.lat[i] <= 60.]
    
    #theta_850.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k')
    theta_equiv_850.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='b')
    #theta_m_equiv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='b', linestyle='-.')
    #theta_m.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    #theta_rc.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    theta_equiv_rc.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='b', linestyle='--')
    theta_equiv_adv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='b', linestyle='-.')
    theta_equiv_net.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C1', linestyle='--')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylim(315.,410.)
    ax.set_xlim(-20.,50.)
    #ax.set_ylim(290.,310.)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('$\Theta$, K')


plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
rcParams['figure.figsize'] = 5.5, 7.
rcParams['font.size'] = 14
    
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)

fig_8('ap_2', ax=ax1, pentad=35)
fig_8('ap_2', ax=ax2, pentad=40)
fig_8('ap_2', ax=ax3, pentad=45)
fig_8('ap_2', ax=ax4, pentad=50)

ax4.set_xlabel('Latitude')

ax1.text(-55, 315., 'a)')
ax2.text(-55, 315., 'b)')
ax3.text(-55, 315., 'c)')
ax4.text(-55, 315., 'd)')

plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)

plt.savefig(plot_dir+'sb08_fig8_moist.pdf', format='pdf')
plt.close()
    


def fig_9(run, ax, pentad=40):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    convTtotheta=(1000./data.pfull)**mc.kappa
    
    heating = data.dt_tg_condensation + data.dt_tg_convection + data.dt_tg_diffusion + data.tdt_rad
    hum_tend = data.dt_qg_condensation + data.dt_qg_convection + data.dt_qg_diffusion
    
    heating_theta_equiv = ((heating + mc.L/mc.cp_air * hum_tend)*convTtotheta).mean('lon')
    
    theta_equiv = (data.temp + mc.L/mc.cp_air * data.sphum/(1-data.sphum)) * convTtotheta
    
    dthetady_equiv = gr.ddy(theta_equiv.mean('lon'), vector=False)
    dthetadp_equiv = gr.ddp(theta_equiv.mean('lon'))
    dthetadt_equiv = gr.ddt(theta_equiv.mean('lon'))
    
    vcomp_theta_equiv = (data.vcomp_temp + mc.L/mc.cp_air * data.sphum_v)*convTtotheta
    vcomp_theta_eddy_equiv = vcomp_theta_equiv.mean('lon') - data.vcomp.mean('lon')*theta_equiv.mean('lon')
    
    vdthetady_mean_equiv = data.vcomp.mean('lon') * dthetady_equiv
    wdthetadp_mean_equiv = data.omega.mean('lon') * dthetadp_equiv
    
    def column_int(var_in):
        var_int = mc.cp_air * var_in.sum('pfull')*5000./mc.grav
        return var_int
    
    vdtdy_mean_int_equiv = -1. * column_int(vdthetady_mean_equiv)
    wdtdp_mean_int_equiv = -1. * column_int(wdthetadp_mean_equiv)
    
    vt_eddy_int_equiv = -1. * column_int(vcomp_theta_eddy_equiv)
    div_vt_eddy_int_equiv = gr.ddy(vt_eddy_int_equiv, vector=True)
    
    heating_theta_int_equiv = column_int(heating_theta_equiv)
    dthetadt_int_equiv = column_int(dthetadt_equiv)
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -60. and data.lat[i] <= 60.]
    
    heating_theta_int_equiv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k')
    dthetadt_int_equiv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='r')
    #wdtdp_mean_int_equiv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C1')
    #vdtdy_mean_int_equiv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C2')
    (vdtdy_mean_int_equiv  + wdtdp_mean_int_equiv ).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    #(-dthetadt_int_equiv + vdtdy_mean_int_equiv + wdtdp_mean_int_equiv + heating_theta_int_equiv + div_vt_eddy_int_equiv).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    (vdtdy_mean_int_equiv + wdtdp_mean_int_equiv + heating_theta_int_equiv + div_vt_eddy_int_equiv).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    div_vt_eddy_int_equiv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylim(-250.,250.)
    #ax.set_xlim(-30.,30.)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('Heating rate, W/m$^2$')
    

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)

fig_9('ap_2', ax=ax1, pentad=35)
fig_9('ap_2', ax=ax2, pentad=40)
fig_9('ap_2', ax=ax3, pentad=45)
fig_9('ap_2', ax=ax4, pentad=50)

ax4.set_xlabel('Latitude')

ax1.text(-55, 500., 'a)')
ax2.text(-55, 500., 'b)')
ax3.text(-55, 500., 'c)')
ax4.text(-55, 500., 'd)')

plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)

plt.savefig(plot_dir+'sb08_fig9_moist.pdf', format='pdf')
plt.close()
