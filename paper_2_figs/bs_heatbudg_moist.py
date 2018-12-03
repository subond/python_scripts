"""
Reproduce moist versions of figures 8 and 9 from SB08 using data from the updated ap_2 run. 
Use Emanuel 1995 (On Thermally Direct Circulations...) angular momentum conserving equivalent 
potential temperature (12/11/2018)
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

    
def get_latmax(data_in):
    data_max = data_in.max('lat')
    data_max_lat = data_in.lat.values[data_in.argmax('lat').values]
    data_max_lat = xr.DataArray(data_max_lat, coords=[data_in.xofyear], dims=['xofyear'])
    return data_max, data_max_lat
        
def fig_8_moist(run, ax, pentad=45, rotfac=1., dayfac=5.):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    heating = data.dt_tg_condensation + data.dt_tg_convection + data.dt_tg_diffusion + data.tdt_rad
    hum_tend = data.dt_qg_condensation + data.dt_qg_convection + data.dt_qg_diffusion
    
    convTtotheta=(1000./data.pfull)**mc.kappa
    theta = data.temp*convTtotheta
    heating_theta = heating*convTtotheta
    heating_equiv_theta = heating_theta + mc.L/mc.cp_air * hum_tend/(1-data.sphum)**2. * convTtotheta
    
    sinphi = np.sin(data.lat * np.pi/180.)
    cosphi = np.cos(data.lat * np.pi/180.)
    
    theta_850 = theta.sel(pfull=850.).mean('lon')
    
    theta_equiv = (data.temp + mc.L/mc.cp_air * data.sphum/(1-data.sphum)) * convTtotheta
    theta_equiv_850 = theta_equiv.sel(pfull=850.).mean('lon')
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -10. and data.lat[i] <= 10.]
    
    Tt = (data.temp.sel(pfull=150.).mean(('lon'))*cosphi).sel(lat=lats).sum('lat')/cosphi.sel(lat=lats).sum('lat')
    Ts = (data.t_surf.mean(('lon'))*cosphi).sel(lat=lats).sum('lat')/cosphi.sel(lat=lats).sum('lat')
    
    dthetady_equiv = gr.ddy(theta_equiv.mean('lon'), vector=False)
    dthetadp_equiv = gr.ddp(theta_equiv.mean('lon'))
    dthetadt_equiv = gr.ddt(theta_equiv.mean('lon'))
    vdthetady_mean = data.vcomp.mean('lon') * dthetady_equiv
    wdthetadp_mean = data.omega.mean('lon') * dthetadp_equiv
    adv_heating_equiv = -1. *(vdthetady_mean + wdthetadp_mean)
    
    w_850 = data.omega.sel(pfull=850.).mean('lon')
    w_max, w_max_lat = get_latmax(-1.*w_850)
    
    theta_max, theta_max_lat = get_latmax(theta_equiv_850)
    cosphimax = np.cos(theta_max_lat * np.pi/180.)
    chi = (rotfac * mc.omega)**2 * mc.a**2./mc.cp_air/(Ts - Tt)
    theta_m_equiv = theta_max * np.exp(-1.* chi * (cosphimax**2.-cosphi**2.)**2./cosphi**2.)

    def adj_theta(theta_in, adj_by, dayfac=dayfac):
        if 'lon' in adj_by.coords:
            return theta_in + adj_by.sel(pfull=850.).mean('lon') * 86400. * dayfac
        else:
            return theta_in + adj_by.sel(pfull=850.) * 86400. * dayfac
    
    theta_equiv_rc = adj_theta(theta_equiv_850, heating_equiv_theta, dayfac=dayfac)
    theta_equiv_adv = adj_theta(theta_equiv_850, adv_heating_equiv, dayfac=dayfac)
    theta_equiv_net = adj_theta(theta_equiv_850, dthetadt_equiv, dayfac=dayfac*3.)
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -60. and data.lat[i] <= 60.]
    
    theta_m_equiv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='0.7')
    theta_equiv_rc.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    theta_equiv_adv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    theta_equiv_850.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C0')
    theta_equiv_net.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C2')
    ax.plot([w_max_lat.sel(xofyear=pentad),w_max_lat.sel(xofyear=pentad)], [300.,380.], color='0.7', linestyle=':')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylim(300.,380.)
    ax.set_xlim(-35.,35.)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('$\Theta$, K')


def plot_bs_fig_8_moist(run, pentads=[35,40,45,50], rotfac=1.):
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/heatbudg_moist/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    rcParams['figure.figsize'] = 5.5, 7.
    rcParams['font.size'] = 14
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    data['precipitation'] = make_sym(data.precipitation)
    precip_centroid(data)       # Locate precipitation centroid
    dpcentdt = gr.ddt(data.p_cent) * 86400.
    p_cent_pos = np.abs(data.p_cent.where(dpcentdt>=0.))
       
    eq_time = p_cent_pos.xofyear[p_cent_pos.argmin('xofyear').values]
    peak_rate_time = dpcentdt.xofyear[dpcentdt.argmax('xofyear').values]
    max_lat_time = data.p_cent.xofyear[data.p_cent.argmax('xofyear').values]
    
    edge_loc, psi_max, psi_max_loc = get_edge_psi(data, lev=850., thresh=0., intdown=True)
    
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
    axes = [ax1,ax2,ax3,ax4]
    pentads = [eq_time.values, peak_rate_time.values, np.floor((peak_rate_time.values + max_lat_time.values)/2.), max_lat_time.values]
    
    for i in range(4):
        edge = edge_loc.sel(xofyear=pentads[i])
        fig_8_moist(run, ax=axes[i], pentad=pentads[i], rotfac=rotfac)
        axes[i].plot([edge, edge], [300., 380.], 'k:')
        axes[i].text(20, 375., 'Pentad ' + str(int(pentads[i])))
    
    ax4.set_xlabel('Latitude')
    
    #ax1.text(20, 315., 'a)')
    #ax2.text(20, 315., 'b)')
    #ax3.text(2, 315., 'c)')
    #ax4.text(-55, 315., 'd)')
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)
    plt.savefig(plot_dir+'sb08_fig8_moist_' + run + '.pdf', format='pdf')
    plt.close()
    


def fig_9_moist(run, ax, pentad=40):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    convTtotheta=(1000./data.pfull)**mc.kappa
    
    heating = data.dt_tg_condensation + data.dt_tg_convection + data.dt_tg_diffusion + data.tdt_rad
    hum_tend = data.dt_qg_condensation + data.dt_qg_convection + data.dt_qg_diffusion
    
    heating_theta_equiv = ((heating + mc.L/mc.cp_air * hum_tend/(1-data.sphum)**2.)*convTtotheta).mean('lon')
    
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
    
    w_850 = data.omega.sel(pfull=850.).mean('lon')
    w_max, w_max_lat = get_latmax(-1.*w_850)
    
    dthetadt_int_equiv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C2')
    heating_theta_int_equiv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    (vdtdy_mean_int_equiv  + wdtdp_mean_int_equiv ).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    #(vdtdy_mean_int_equiv + wdtdp_mean_int_equiv + heating_theta_int_equiv + div_vt_eddy_int_equiv).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    div_vt_eddy_int_equiv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    ax.plot([w_max_lat.sel(xofyear=pentad),w_max_lat.sel(xofyear=pentad)], [-250.,250.], color='0.7', linestyle=':')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylim(-250.,250.)
    ax.set_xlim(-35.,35.)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('Heating rate, W/m$^2$')
    


def plot_bs_fig_9_moist(run, pentads=[35,40,45,50], rotfac=1.):
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/heatbudg_moist/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    rcParams['figure.figsize'] = 5.5, 7.
    rcParams['font.size'] = 14
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    data['precipitation'] = make_sym(data.precipitation)
    precip_centroid(data)       # Locate precipitation centroid
    dpcentdt = gr.ddt(data.p_cent) * 86400.
    p_cent_pos = np.abs(data.p_cent.where(dpcentdt>=0.))
       
    eq_time = p_cent_pos.xofyear[p_cent_pos.argmin('xofyear').values]
    peak_rate_time = dpcentdt.xofyear[dpcentdt.argmax('xofyear').values]
    max_lat_time = data.p_cent.xofyear[data.p_cent.argmax('xofyear').values]
    
    edge_loc, psi_max, psi_max_loc = get_edge_psi(data, lev=850., thresh=0., intdown=True)
    
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
    axes = [ax1,ax2,ax3,ax4]
    pentads = [eq_time.values, peak_rate_time.values, np.floor((peak_rate_time.values + max_lat_time.values)/2.), max_lat_time.values]
    
    for i in range(4):
        fig_9_moist(run, ax=axes[i], pentad=pentads[i])
        edge = edge_loc.sel(xofyear=pentads[i])
        axes[i].plot([edge, edge], [-250., 250.], 'k:')        
        axes[i].text(20, 240., 'Pentad ' + str(int(pentads[i])))
    
    ax4.set_xlabel('Latitude')
    
    #ax1.text(-55, 315., 'a)')
    #ax2.text(-55, 315., 'b)')
    #ax3.text(-55, 315., 'c)')
    #ax4.text(-55, 315., 'd)')
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)
    plt.savefig(plot_dir+'sb08_fig9_moist_' + run + '.pdf', format='pdf')
    plt.close()
    

plot_bs_fig_8_moist('mld_2.5_heatbudg')
plot_bs_fig_8_moist('mld_5_heatbudg')
plot_bs_fig_8_moist('mld_15_heatbudg')
plot_bs_fig_8_moist('mld_20_heatbudg')

plot_bs_fig_9_moist('mld_2.5_heatbudg')
plot_bs_fig_9_moist('mld_5_heatbudg')
plot_bs_fig_9_moist('mld_15_heatbudg')
plot_bs_fig_9_moist('mld_20_heatbudg')
    
#plot_bs_fig_8_moist('sn_1.000_evap_fluxes_heattrans')
#plot_bs_fig_8_moist('rt_2.000_heatbudg', rotfac=2.)
#plot_bs_fig_8_moist('rt_0.500_heatbudg', rotfac=0.5)
#plot_bs_fig_8_moist('rt_0.750_heatbudg', rotfac=0.75)
#plot_bs_fig_8_moist('rt_1.250_heatbudg', rotfac=1.25)
#plot_bs_fig_8_moist('rt_1.500_heatbudg', rotfac=1.5)
#plot_bs_fig_8_moist('rt_1.750_heatbudg', rotfac=1.75)

#plot_bs_fig_9_moist('sn_1.000_evap_fluxes_heattrans')
#plot_bs_fig_9_moist('rt_2.000_heatbudg')
#plot_bs_fig_9_moist('rt_0.500_heatbudg')
#plot_bs_fig_9_moist('rt_0.750_heatbudg')
#plot_bs_fig_9_moist('rt_1.250_heatbudg')
#plot_bs_fig_9_moist('rt_1.500_heatbudg')
#plot_bs_fig_9_moist('rt_1.750_heatbudg')