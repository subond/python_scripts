"""
Create momentum budget version of figure 8 from SB08 using data from the updated sn_1.000 run (01/09/2018)
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


def bs_mombudg(run, ax, pentad=45, lev=150., dayfac=3., rotfac=1.):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    damping = data.dt_ug_diffusion# + data.tdt_diss_rdamp  
    
    dudy = gr.ddy(data.ucomp.mean('lon'))
    duvdy = gr.ddy(data.ucomp_vcomp.mean('lon'), uv=True)
    dudp = gr.ddp(data.ucomp.mean('lon'))
    duwdp = gr.ddp(data.ucomp_omega.mean('lon'))
    dudt = gr.ddt(data.ucomp.mean('lon'))
    
    vdudy_mean = data.vcomp.mean('lon') * dudy
    wdudp_mean = data.omega.mean('lon') * dudp
    adv_tend = -1. *(vdudy_mean + wdudp_mean)
    vert_adv = -1.*wdudp_mean
    horiz_adv = -1.*vdudy_mean
    
    eddy_tend = -1.*duvdy -1.*duwdp - adv_tend
    vert_eddy = -1.*duwdp - vert_adv
    horiz_eddy = -1.*duvdy - horiz_adv
    
    sinphi = np.sin(data.lat * np.pi/180.)
    cosphi = np.cos(data.lat * np.pi/180.)
    
    coriolis = 2.* rotfac * mc.omega * sinphi * data.vcomp
    
    
    def get_latmax(data_in):
        data_max = data_in.max('lat')
        data_max_lat = data_in.lat.values[data_in.argmax('lat').values]
        data_max_lat = xr.DataArray(data_max_lat, coords=[data_in.xofyear], dims=['xofyear'])
        return data_max, data_max_lat
    
    w_850 = data.omega.sel(pfull=850.).mean('lon')
    w_max, w_max_lat = get_latmax(-1.*w_850)
    convTtotheta=(1000./data.pfull)**mc.kappa
    theta_equiv_850 = ((data.temp + mc.L/mc.cp_air * data.sphum/(1-data.sphum)) * convTtotheta).sel(pfull=850.).mean('lon')
    theta_max, theta_max_lat = get_latmax(theta_equiv_850)
    
    cosphimax = np.cos(theta_max_lat * np.pi/180.)
    
    u_150 = data.ucomp.sel(pfull=lev).mean('lon')
    u_m = rotfac * mc.omega * mc.a * (cosphimax**2. - cosphi**2.)/cosphi
    
    def adj_u(u_in, adj_by, dayfac=dayfac):
        if 'lon' in adj_by.coords:
            return u_in + adj_by.sel(pfull=lev).mean('lon') * 86400. * dayfac
        else:
            return u_in + adj_by.sel(pfull=lev) * 86400. * dayfac
    
    u_damp = adj_u(u_150, damping)
    u_cor = adj_u(u_150, coriolis)
    u_adv = adj_u(u_150, adv_tend);   u_vadv = adj_u(u_150, horiz_adv);    u_hadv = adj_u(u_150, vert_adv)
    u_eddy = adj_u(u_150, eddy_tend);   u_veddy = adj_u(u_150, horiz_eddy);    u_heddy = adj_u(u_150, vert_eddy)
    
    resid = coriolis + damping + adv_tend + eddy_tend
    u_resid = adj_u(u_150, resid, dayfac=10.*dayfac)
    
    u_net = adj_u(u_150, dudt, dayfac=10.*dayfac)
    
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -60. and data.lat[i] <= 60.]
    
    u_m.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='0.7')
    u_cor.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='--')
    u_adv.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle='-.')
    u_eddy.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='k', linestyle=':')
    
    u_150.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C0')
    u_net.sel(xofyear=pentad, lat=lats).plot(ax=ax, color='C2')
    ax.plot([w_max_lat.sel(xofyear=pentad),w_max_lat.sel(xofyear=pentad)], [-75.,75.], color='0.7', linestyle=':')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_xlim(-35.,35.)
    ax.set_ylim(-75.,75.)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('u, m/s')



def plot_bs_mombudg(run, pentads=[35,40,45,50], rotfac=1.):
    plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/mombudg/'
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
        bs_mombudg(run, ax=axes[i], pentad=pentads[i], rotfac=rotfac)
        axes[i].text(20, 70., 'Pentad ' + str(int(pentads[i])))
        edge = edge_loc.sel(xofyear=pentads[i])
        axes[i].plot([edge, edge], [-75., 75.], 'k:')
    
    ax4.set_xlabel('Latitude')
    #ax1.text(-55, 315., 'a)')
    #ax2.text(-55, 315., 'b)')
    #ax3.text(-55, 315., 'c)')
    #ax4.text(-55, 315., 'd)')
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.1)
    plt.savefig(plot_dir+'sb08_mombudg_' + run + '.pdf', format='pdf')
    plt.close()


plot_bs_mombudg('mld_2.5_heatbudg')
plot_bs_mombudg('mld_5_heatbudg')
plot_bs_mombudg('mld_15_heatbudg')
plot_bs_mombudg('mld_20_heatbudg')

#plot_bs_mombudg('sn_1.000_evap_fluxes_heattrans')
#plot_bs_mombudg('rt_2.000_heatbudg', rotfac=2.)
#plot_bs_mombudg('rt_0.500_heatbudg', rotfac=0.5)
#plot_bs_mombudg('rt_0.750_heatbudg', rotfac=0.75)
#plot_bs_mombudg('rt_1.250_heatbudg', rotfac=1.25)
#plot_bs_mombudg('rt_1.500_heatbudg', rotfac=1.5)
#plot_bs_mombudg('rt_1.750_heatbudg', rotfac=1.75)