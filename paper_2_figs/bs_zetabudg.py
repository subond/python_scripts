"""
Create vorticity budget version of figure 8 from SB08 using data from the updated sn_1.000 run (01/09/2018)
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


def fig_8(run, ax, pentad=45, lev=200., dayfac=3., rotfac=1.):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    sinphi = np.sin(data.lat * np.pi/180.)
    cosphi = np.cos(data.lat * np.pi/180.)
        
    zeta = 2.* rotfac * mc.omega*sinphi -1.* gr.ddy(data.ucomp.mean('lon')) #* 86400.
    nu = gr.ddp(data.ucomp.mean('lon')) #* 86400.
    
    dzetady = gr.ddy(zeta, vector=False)
    dzetadp = gr.ddp(zeta)
    dwdy = gr.ddy(data.omega, vector=False)
    
    dvdy = gr.ddy(data.vcomp)
    dwdp = gr.ddp(data.omega)
    
    dzetadt = gr.ddt(zeta)
    
    vdzetady_mean = data.vcomp.mean('lon') * dzetady
    wdzetadp_mean = data.omega.mean('lon') * dzetadp
    adv_tend = -1. *(vdzetady_mean + wdzetadp_mean)
    vert_adv = -1.*wdzetadp_mean
    horiz_adv = -1.*vdzetady_mean
    
    stretching = -1.* zeta * dvdy
    tilting = nu * dwdy
    
    eddy_tend = dzetadt - adv_tend - stretching - tilting
    
        
    w_850 = data.omega.sel(pfull=850.).mean('lon')
    w_max = w_850.min('lat')
    w_max_lat = data.lat.values[w_850.argmin('lat').values]
    w_max_lat = xr.DataArray(w_max_lat, coords=[data.xofyear], dims=['xofyear'])
    cosphimax = np.cos(w_max_lat * np.pi/180.)
    
    u_m = rotfac * mc.omega * mc.a * (cosphimax**2. - cosphi**2.)/cosphi
    zeta_m = 2.*mc.omega*sinphi -1.* gr.ddy(u_m)
    
    def adj_zeta(zeta_in, adj_by, dayfac=dayfac):
        if 'lon' in adj_by.coords:
            return zeta_in + adj_by.sel(pfull=lev).mean('lon') * 86400.**2. * dayfac
        else:
            return zeta_in + adj_by.sel(pfull=lev) * 86400.**2. * dayfac
    
    zeta_150 = zeta.sel(pfull=150.)*86400.
    
    zeta_adv = adj_zeta(zeta_150, adv_tend);   zeta_vadv = adj_zeta(zeta_150, horiz_adv);    zeta_hadv = adj_zeta(zeta_150, vert_adv)
    zeta_stretching = adj_zeta(zeta_150, stretching)
    zeta_tilting = adj_zeta(zeta_150, tilting)
    zeta_eddy = adj_zeta(zeta_150, eddy_tend);
    
    zeta_net = adj_zeta(zeta_150, dzetadt, dayfac=10.*dayfac)
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -60. and data.lat[i] <= 60.]
    
    zeta_m.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='0.7', linestyle='-')
    zeta_stretching.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    zeta_tilting.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='0.7', linestyle='--')

    zeta_adv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    #u_hadv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C1', linestyle='-.')
    #u_vadv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C2', linestyle='-.')
    zeta_eddy.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    
    zeta_150.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C0')
    zeta_net.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C2')
    ax.plot([w_max_lat.sel(xofyear=pentad),w_max_lat.sel(xofyear=pentad)], [-7.5,7.5], color='0.7', linestyle=':')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_xlim(-35.,35.)
    ax.set_ylim(-7.5,7.5)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('zeta, /s')


def plot_bs_zetabudg(run, pentads=[35,40,45,55], rotfac=1.):
    plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/zetabudg/'
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
        axes[i].text(20, 7.0, 'Pentad ' + str(int(pentads[i])))
        edge = edge_loc.sel(xofyear=pentads[i])
        axes[i].plot([edge, edge], [-7.5, 7.5], 'k:')
    
    ax4.set_xlabel('Latitude')
    
    #ax1.text(-55, 315., 'a)')
    #ax2.text(-55, 315., 'b)')
    #ax3.text(-55, 315., 'c)')
    #ax4.text(-55, 315., 'd)')
    plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)
    plt.savefig(plot_dir+'sb08_zetabudg_fig8_' + run + '.pdf', format='pdf')
    plt.close()

plot_bs_zetabudg('mld_2.5_heatbudg')
plot_bs_zetabudg('mld_5_heatbudg')
plot_bs_zetabudg('mld_15_heatbudg')
plot_bs_zetabudg('mld_20_heatbudg')


#plot_bs_zetabudg('sn_1.000_evap_fluxes_heattrans')
#plot_bs_zetabudg('rt_2.000_heatbudg', rotfac=2.)
#plot_bs_zetabudg('rt_0.500_heatbudg', rotfac=0.5)
#plot_bs_zetabudg('rt_0.750_heatbudg', rotfac=0.75)
#plot_bs_zetabudg('rt_1.250_heatbudg', rotfac=1.25)
#plot_bs_zetabudg('rt_1.500_heatbudg', rotfac=1.5)
#plot_bs_zetabudg('rt_1.750_heatbudg', rotfac=1.75)
#plot_bs_zetabudg('sn_1_sst_zs')
