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
    
    zeta = 2.*rotfac * mc.omega*sinphi -1.* gr.ddy(data.ucomp.mean('lon')) #* 86400.
    nu = gr.ddp(data.ucomp.mean('lon')) #* 86400.
    
    dnudy = gr.ddy(nu, vector=False)
    dnudp = gr.ddp(nu)
    dwdp = gr.ddp(data.omega)
    dvdp = gr.ddp(data.vcomp)
    
    dnudt = gr.ddt(nu)
    
    vdnudy_mean = data.vcomp.mean('lon') * dnudy
    wdnudp_mean = data.omega.mean('lon') * dnudp
    adv_tend = -1. *(vdnudy_mean + wdnudp_mean)
    vert_adv = -1.*wdnudp_mean
    horiz_adv = -1.*vdnudy_mean
    
    stretching = -1.* nu * dwdp
    tilting = zeta * dvdp
    
    eddy_tend = dnudt - adv_tend - stretching - tilting
    
        
    w_850 = data.omega.sel(pfull=850.).mean('lon')
    w_max = w_850.min('lat')
    w_max_lat = data.lat.values[w_850.argmin('lat').values]
    w_max_lat = xr.DataArray(w_max_lat, coords=[data.xofyear], dims=['xofyear'])
    cosphimax = np.cos(w_max_lat * np.pi/180.)
    
    #u_m = mc.omega * mc.a * (cosphimax**2. - cosphi**2.)/cosphi
    #nu_m = 2.*mc.omega*sinphi -1.* gr.ddy(u_m)
    
    def adj_nu(nu_in, adj_by, dayfac=dayfac):
        if 'lon' in adj_by.coords:
            return nu_in + adj_by.sel(pfull=lev).mean('lon') * 86400.**2. * dayfac
        else:
            return nu_in + adj_by.sel(pfull=lev) * 86400.**2. * dayfac
    
    nu_150 = nu.sel(pfull=150.)*86400.
    
    nu_adv = adj_nu(nu_150, adv_tend);   nu_vadv = adj_nu(nu_150, horiz_adv);    nu_hadv = adj_nu(nu_150, vert_adv)
    nu_stretching = adj_nu(nu_150, stretching)
    nu_tilting = adj_nu(nu_150, tilting)
    nu_eddy = adj_nu(nu_150, eddy_tend);
    
    nu_net = adj_nu(nu_150, dnudt, dayfac=10.*dayfac)
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -30. and data.lat[i] <= 30.]
    
    #nu_m.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='0.7', linestyle='-')
    nu_stretching.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    nu_tilting.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='0.7', linestyle='--')

    nu_adv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    #u_hadv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C1', linestyle='-.')
    #u_vadv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C2', linestyle='-.')
    nu_eddy.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    
    nu_150.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C0')
    nu_net.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C2')
    ax.plot([w_max_lat.sel(xofyear=pentad),w_max_lat.sel(xofyear=pentad)], [-300,300], color='0.7', linestyle=':')
    ax.set_title('')
    ax.set_xlabel('')
    #ax.set_ylim(285.,315.)
    ax.set_ylim(-300,300)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('nu, /s')


def plot_bs_nubudg(run, pentads=[35,40,45,55], rotfac=1.):
    plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/nubudg/'
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
        fig_8(run, ax=axes[i], pentad=pentads[i], rotfac=1.)
        axes[i].text(20, 290., 'Pentad ' + str(int(pentads[i])))
        edge = edge_loc.sel(xofyear=pentads[i])
        axes[i].plot([edge, edge], [-300., 300.], 'k:')
    
    ax4.set_xlabel('Latitude')
    
    #ax1.text(-55, 315., 'a)')
    #ax2.text(-55, 315., 'b)')
    #ax3.text(-55, 315., 'c)')
    #ax4.text(-55, 315., 'd)')
    plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)
    plt.savefig(plot_dir+'sb08_nubudg_fig8_' + run + '.pdf', format='pdf')
    plt.close()

plot_bs_nubudg('mld_2.5_heatbudg')
plot_bs_nubudg('mld_5_heatbudg')
plot_bs_nubudg('mld_15_heatbudg')
plot_bs_nubudg('mld_20_heatbudg')

#plot_bs_nubudg('sn_1.000_evap_fluxes_heattrans')
#plot_bs_nubudg('rt_2.000_heatbudg', rotfac=2.)
#plot_bs_nubudg('rt_0.500_heatbudg', rotfac=0.5)
#plot_bs_nubudg('rt_0.750_heatbudg', rotfac=0.75)
#plot_bs_nubudg('rt_1.250_heatbudg', rotfac=1.25)
#plot_bs_nubudg('rt_1.500_heatbudg', rotfac=1.5)
#plot_bs_nubudg('rt_1.750_heatbudg', rotfac=1.75)
#plot_bs_nubudg('sn_1_sst_zs')
