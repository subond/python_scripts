"""
Reproduce figure 4 from SB08 using data for mld and rot runs (21/11/2018)
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

        
def fig_4(run, rotfac=1.):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    psi_full = mass_streamfunction(data, a=6376.0e3, dp_in=50.) / 1.e9
    vcomp = data.vcomp
    
    sinphi = np.sin(data.lat * np.pi/180.)
    cosphi = np.cos(data.lat * np.pi/180.)
    coriolis = 2.* rotfac * mc.omega * sinphi 
    
    dudy = gr.ddy(data.ucomp)
    duvdy = gr.ddy(data.ucomp_vcomp, uv=True)
    dudp = gr.ddp(data.ucomp)
    duwdp = gr.ddp(data.ucomp_omega)
    
    vdudy_mean = data.vcomp * dudy
    wdudp_mean = data.omega * dudp
    adv_v = (vdudy_mean + wdudp_mean)/coriolis
    eddy_v = (duvdy + duwdp)/coriolis - adv_v
    
    def psi_comp(comp, v):
        data['vcomp'] = comp
        psi_comp = mass_streamfunction(data, a=6376.0e3, dp_in=50.) / 1.e9
        data['vcomp'] = v
        return psi_comp
    
    psi_eddy = psi_comp(eddy_v, vcomp)
    psi_adv = psi_comp(adv_v, vcomp)
    
    
    # Get times to plot at - time when on equator, when centroid is moving fastest, when off the equator
    data['precipitation'] = make_sym(data.precipitation)
    precip_centroid(data)       # Locate precipitation centroid
    dpcentdt = gr.ddt(data.p_cent) * 86400.
    p_cent_pos = np.abs(data.p_cent.where(dpcentdt>=0.))
       
    eq_time = p_cent_pos.xofyear[p_cent_pos.argmin('xofyear').values]
    peak_rate_time = dpcentdt.xofyear[dpcentdt.argmax('xofyear').values]
    max_lat_time = data.p_cent.xofyear[data.p_cent.argmax('xofyear').values]
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/cell_components/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    rcParams['figure.figsize'] = 7., 7.  # X, Y
    rcParams['font.size'] = 14
        
    fig, ((ax1,ax2,ax3), (ax4,ax5,ax6), (ax7,ax8,ax9), (ax10,ax11,ax12)) = plt.subplots(4, 3, sharex=True, sharey=True)
    axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12]
    pentads = [eq_time.values, peak_rate_time.values, np.floor((peak_rate_time.values + max_lat_time.values)/2.), max_lat_time.values]
    
    levels = np.arange(-500.,501.,50.)
    levels_pos = np.arange(0.,601.,50.)
    levels_neg = np.arange(-600.,0.,50.)
    for i in range(4):
        psi_full.sel(xofyear=pentads[i]).plot.contourf(x='lat', y='pfull', yincrease=False, ax=axes[i*3], levels=levels, add_labels=False, add_colorbar=False)
        #psi_full.sel(xofyear=pentads[i]).plot.contour(x='lat', y='pfull', yincrease=False, ax=axes[i*3], levels=levels_pos, add_labels=False, colors='k', linestyles='--')
        #psi_full.sel(xofyear=pentads[i]).plot.contour(x='lat', y='pfull', yincrease=False, ax=axes[i*3], levels=levels_neg, add_labels=False, colors='k')
        psi_adv.sel(xofyear=pentads[i]).plot.contourf(x='lat', y='pfull', yincrease=False, ax=axes[i*3+1], levels=levels, add_labels=False, add_colorbar=False)
        #psi_adv.sel(xofyear=pentads[i]).plot.contour(x='lat', y='pfull', yincrease=False, ax=axes[i*3+1], levels=levels_pos, add_labels=False, colors='k', linestyles='--')
        #psi_adv.sel(xofyear=pentads[i]).plot.contour(x='lat', y='pfull', yincrease=False, ax=axes[i*3+1], levels=levels_neg, add_labels=False, colors='k')
        psi_eddy.sel(xofyear=pentads[i]).plot.contourf(x='lat', y='pfull', yincrease=False, ax=axes[i*3+2], levels=levels, add_labels=False, add_colorbar=False)
        #psi_eddy.sel(xofyear=pentads[i]).plot.contour(x='lat', y='pfull', yincrease=False, ax=axes[i*3+2], levels=levels_pos, add_labels=False, colors='k', linestyles='--')
        #psi_eddy.sel(xofyear=pentads[i]).plot.contour(x='lat', y='pfull', yincrease=False, ax=axes[i*3+2], levels=levels_neg, add_labels=False, colors='k')
    
    for ax in axes:
        ax.fill_between([-2.5,2.5], 0., 850.,  facecolor='white')
        ax.set_xlim(-40.,40.)
        ax.set_ylim(850.,0.)
        ax.set_yticks([800.,600.,400.,200.])
        #ax.text(20, 310., 'Pentad ' + str(int(pentads[i])))
    
    ax1.set_title('Full $\Psi$')
    ax2.set_title('Mean state $\Psi$')
    ax3.set_title('Eddy $\Psi$')
    ax10.set_xlabel('Latitude')
    ax11.set_xlabel('Latitude')
    ax12.set_xlabel('Latitude')
    ax1.set_ylabel('Pressure, hPa')
    ax4.set_ylabel('Pressure, hPa')
    ax7.set_ylabel('Pressure, hPa')
    ax10.set_ylabel('Pressure, hPa')
    
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
    plt.savefig(plot_dir+'sb08_fig4_' + run + '.pdf', format='pdf')
    plt.close()


fig_4('sn_1.000_evap_fluxes_heattrans')
fig_4('rt_0.500', rotfac=0.5)
fig_4('rt_0.750', rotfac=0.75)
fig_4('rt_1.250', rotfac=1.25)
fig_4('rt_1.500', rotfac=1.5)
fig_4('rt_1.750', rotfac=1.75)
fig_4('rt_2.000', rotfac=2.)
fig_4('mld_2.5')
fig_4('mld_5')
fig_4('mld_15')
fig_4('mld_20')
