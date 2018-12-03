"""
Reproduce figures 8 and 9 from SB08 using data from the updated sn_1.000 run (01/09/2018)
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

    
def get_latmax(data_in):
    data_max = data_in.max('lat')
    data_max_lat = data_in.lat.values[data_in.argmax('lat').values]
    data_max_lat = xr.DataArray(data_max_lat, coords=[data_in.xofyear], dims=['xofyear'])
    return data_max, data_max_lat
        
def fig_8(run, ax, pentad=45, rotfac=1., dayfac=1.):

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
    
    
    w_850 = data.omega.sel(pfull=850.).mean('lon')
    w_max, w_max_lat = get_latmax(-1.*w_850)
    theta_850 = theta.sel(pfull=850.).mean('lon')
    theta_max, theta_max_lat = get_latmax(theta_850)
    
    sinphimax = np.sin(theta_max_lat * np.pi/180.)
    
    theta_m = theta_max - 300.*(mc.omega * rotfac)**2.*mc.a**2./(2.*mc.grav*14000.) * (sinphi**2.-sinphimax**2.)**2./cosphi**2.
    
    def adj_theta(theta_in, adj_by, dayfac=dayfac):
        if 'lon' in adj_by.coords:
            return theta_in + adj_by.sel(pfull=850.).mean('lon') * 86400. * dayfac
        else:
            return theta_in + adj_by.sel(pfull=850.) * 86400. * dayfac
    
    
    theta_rc = adj_theta(theta_850, heating_theta)
    theta_adv = adj_theta(theta_850, adv_heating)
    theta_net = adj_theta(theta_850, dthetadt, dayfac=20.)
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -60. and data.lat[i] <= 60.]
    
    theta_m.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='0.7')
    theta_rc.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    theta_adv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    theta_850.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C0')
    theta_net.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C2')
    ax.plot([w_max_lat.sel(xofyear=pentad),w_max_lat.sel(xofyear=pentad)], [295,315], color='0.7', linestyle=':')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylim(295.,315.)
    ax.set_xlim(-35.,35.)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('$\Theta$, K')


def plot_bs_fig_8(run, pentads=[35,40,45,50], rotfac=1.):
    plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/heatbudg/'
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
        fig_8(run, ax=axes[i], pentad=pentads[i], rotfac=rotfac)
        axes[i].text(20, 310., 'Pentad ' + str(int(pentads[i])))
        edge = edge_loc.sel(xofyear=pentads[i])
        axes[i].plot([edge, edge], [295., 315.], 'k:')
    
    ax4.set_xlabel('Latitude')
    
    #ax1.text(-55, 315., 'a)')
    #ax2.text(-55, 315., 'b)')
    #ax3.text(-55, 315., 'c)')
    #ax4.text(-55, 315., 'd)')
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)
    plt.savefig(plot_dir+'sb08_fig8_' + run + '.pdf', format='pdf')
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
    
    
    w_850 = data.omega.sel(pfull=850.).mean('lon')
    w_max, w_max_lat = get_latmax(-1.*w_850)
    
    heating_theta_int = column_int(heating_theta)
    dthetadt_int = column_int(dthetadt)
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -40. and data.lat[i] <= 40.]
    
    dthetadt_int.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C2')
    heating_theta_int.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    #wdtdp_mean_int.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C1')
    #vdtdy_mean_int.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    (vdtdy_mean_int + wdtdp_mean_int).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    #(-dthetadt_int + vdtdy_mean_int + wdtdp_mean_int + heating_theta_int + div_vt_eddy_int).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    #(vdtdy_mean_int + wdtdp_mean_int + heating_theta_int + div_vt_eddy_int).sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    div_vt_eddy_int.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    ax.plot([w_max_lat.sel(xofyear=pentad),w_max_lat.sel(xofyear=pentad)], [-500.,500.], color='0.7', linestyle=':')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylim(-500.,500.)
    ax.set_xlim(-35.,35.)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('Heating rate, W/m$^2$')
    

def plot_bs_fig_9(run, pentads=[35,40,45,50]):
    plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/heatbudg/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
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
        fig_9(run, ax=axes[i], pentad=pentads[i])
        axes[i].text(20, 315., 'Pentad ' + str(int(pentads[i])))
        edge = edge_loc.sel(xofyear=pentads[i])
        axes[i].plot([edge, edge], [-500., 500.], 'k:')
    
    ax4.set_xlabel('Latitude')
    
    ax1.text(-55, 500., 'a)')
    ax2.text(-55, 500., 'b)')
    ax3.text(-55, 500., 'c)')
    ax4.text(-55, 500., 'd)')
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)
    plt.savefig(plot_dir+'sb08_fig9_' + run + '.pdf', format='pdf')
    plt.close()


plot_bs_fig_8('mld_2.5_heatbudg')
plot_bs_fig_8('mld_5_heatbudg')
plot_bs_fig_8('mld_15_heatbudg')
plot_bs_fig_8('mld_20_heatbudg')

plot_bs_fig_9('mld_2.5_heatbudg')
plot_bs_fig_9('mld_5_heatbudg')
plot_bs_fig_9('mld_15_heatbudg')
plot_bs_fig_9('mld_20_heatbudg')

#plot_bs_fig_8('sn_1.000_evap_fluxes_heattrans')
#plot_bs_fig_8('rt_2.000_heatbudg', rotfac=2.)
#plot_bs_fig_8('rt_0.500_heatbudg', rotfac=0.5)
#plot_bs_fig_8('rt_0.750_heatbudg', rotfac=0.75)
#plot_bs_fig_8('rt_1.250_heatbudg', rotfac=1.25)
#plot_bs_fig_8('rt_1.500_heatbudg', rotfac=1.5)
#plot_bs_fig_8('rt_1.750_heatbudg', rotfac=1.75)

#plot_bs_fig_9('sn_1.000_evap_fluxes_heattrans')
#plot_bs_fig_9('rt_2.000_heatbudg')
#plot_bs_fig_9('rt_0.500_heatbudg')
#plot_bs_fig_9('rt_0.750_heatbudg')
#plot_bs_fig_9('rt_1.250_heatbudg')
#plot_bs_fig_9('rt_1.500_heatbudg')
#plot_bs_fig_9('rt_1.750_heatbudg')
