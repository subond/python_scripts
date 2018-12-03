"""
Code up the elliptical equation for the mass streamfunction and plot for when the ITCZ is on the equator, the fastest moving pentad, and the furthest from the equator (20/11/2018)
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


def elliptical_equation(run, rotfac=1.):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    data = data.mean('lon')
    
    data['precipitation'] = make_sym(data.precipitation)
    precip_centroid(data)       # Locate precipitation centroid
    dpcentdt = gr.ddt(data.p_cent) * 86400.
    p_cent_pos = np.abs(data.p_cent.where(dpcentdt>=0.))
       
    eq_time = p_cent_pos.xofyear[p_cent_pos.argmin('xofyear').values]
    peak_rate_time = dpcentdt.xofyear[dpcentdt.argmax('xofyear').values]
    max_lat_time = data.p_cent.xofyear[data.p_cent.argmax('xofyear').values]
    
    convTtotheta=(1000./data.pfull)**mc.kappa
    data['theta'] = data.temp*convTtotheta 
    #data['vcomp_theta'] = data.vcomp_temp*convTtotheta 
    #data['omega_theta'] = data.omega_temp*convTtotheta 
    #heating = data.dt_tg_condensation + data.dt_tg_convection + data.dt_tg_diffusion + data.tdt_rad
    #data['heating_theta'] = heating*convTtotheta
    
    sinphi = np.sin(data.lat * np.pi/180.)
    cosphi = np.cos(data.lat * np.pi/180.)
    
    #Calculate psi 
    #psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    
    dpsidydy = gr.ddy(data.omega, vector=False)
    rho = data.pfull*100./mc.rdgas/data.temp
    dthetadp = gr.ddp(data.theta)
    brunt_fac = -1./rho/data.theta * dthetadp
    psi_yy = (brunt_fac * dpsidydy)
    
    dpsidpdp = gr.ddp(-1.*data.vcomp)
    coriolis = 2.* rotfac * mc.omega * sinphi * data.vcomp
    dudy = gr.ddy(data.ucomp)
    cor_fac = -1.*dudy*coriolis + coriolis**2.
    psi_pp = (cor_fac * dpsidydy)
    
    dpsidpdy = gr.ddy(-1.*data.vcomp)
    dudp = gr.ddp(data.ucomp)
    psi_py = (2.*dudp*dpsidpdy*coriolis)
    
    duvdy = gr.ddy(data.ucomp_vcomp, uv=True)
    duwdp = gr.ddp(data.ucomp_omega)
    vdudy_mean = data.vcomp * dudy
    wdudp_mean = data.omega * dudp
    eddy_tend = -1.*duvdy -1.*duwdp + vdudy_mean + wdudp_mean
    M = data.dt_ug_diffusion + eddy_tend
    dMdp = (gr.ddp(M) * coriolis)
    
    dthetady = gr.ddy(data.theta, vector=False)
    #dvtdy = gr.ddy(data.vcomp_theta)
    #dwtdp = gr.ddp(data.omega_theta)
    vdtdy_mean = data.vcomp * dthetady
    wdtdp_mean = data.omega * dthetadp
    #eddy_tend_th = -1.*dvtdy -1.*dwtdp + vdtdy_mean + wdtdp_mean
    #J = data.heating_theta #+ eddy_tend_th
    #dJdy = (-1.* gr.ddy(J, vector=False)/rho/data.theta)
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/elliptical_eq/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    rcParams['figure.figsize'] = 20., 12.
    
    levels = np.arange(-1.e-12,1.1e-12,0.1e-12)
    fig, ((ax1, ax2, ax3, ax4, ax5), (ax6, ax7, ax8, ax9, ax10), (ax11, ax12, ax13, ax14, ax15)) = plt.subplots(3, 5, sharey='row', sharex='col')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15]
    
    i=0
    for time in [eq_time, peak_rate_time, max_lat_time]:
        f1 = psi_yy.sel(xofyear=time).plot.contourf(ax=axes[5*i], x='lat', y='pfull', yincrease=False, levels=levels, add_labels=False, extend='both', add_colorbar=False)
        psi_pp.sel(xofyear=time).plot.contourf(ax=axes[5*i+1], x='lat', y='pfull', yincrease=False, levels=levels, add_labels=False, extend='both', add_colorbar=False)
        psi_py.sel(xofyear=time).plot.contourf(ax=axes[5*i+2], x='lat', y='pfull', yincrease=False, levels=levels, add_labels=False, extend='both', add_colorbar=False)
        dMdp.sel(xofyear=time).plot.contourf(ax=axes[5*i+3], x='lat', y='pfull', yincrease=False, levels=levels, add_labels=False, extend='both', add_colorbar=False)
        #dJdy.sel(xofyear=time).plot.contourf(ax=axes[5*i+4], x='lat', y='pfull', yincrease=False, levels=levels, add_labels=False, extend='both', add_colorbar=False)
        i=i+1
    for ax in axes:
        ax.set_xlim(-30.,30.)
    
    ax1.set_title('N$^2d^2psi/dy^2$')
    ax2.set_title('$f(f-du/dy)d^2psi/dp^2$')
    ax3.set_title('Cross term')
    ax4.set_title('dM/dp')
    ax5.set_title('dJ/dy')
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.25, wspace=0.12)
    #Colorbar
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    
    figname = 'elliptical_eq_' +run+ '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()    

#elliptical_equation('sn_1.000_evap_fluxes_heattrans')
elliptical_equation('sn_1_sst_zs')
