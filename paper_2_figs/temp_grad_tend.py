"""
Make lat-time plots of the terms modifying the meridional potential temperature/equivalent potential temperature gradient
Neglect eddies to start (16/11/2018)
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

    
def temp_grad_tend(run, lev=850., moist=False):
    
    rcParams['figure.figsize'] = 12, 5.5
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    psi = np.abs(psi).max('pfull')
    
    data = data.mean('lon')
    convTtotheta=(1000./data.pfull)**mc.kappa
    heating = (data.dt_tg_condensation + data.dt_tg_convection + data.dt_tg_diffusion + data.tdt_rad)
    hum_tend = data.dt_qg_condensation + data.dt_qg_convection + data.dt_qg_diffusion
    
    if moist:
        heating = (heating + mc.L/mc.cp_air * hum_tend/(1-data.sphum)**2.) * convTtotheta    
        theta = (data.temp + mc.L/mc.cp_air * data.sphum/(1-data.sphum)) * convTtotheta
    else:
        theta = data.temp*convTtotheta
        #vcomp_theta = data.vcomp_temp*convTtotheta
        #omega_theta = data.omega_temp*convTtotheta
        heating = heating*convTtotheta
    
    dthetady = gr.ddy(theta, vector=False)
    dthetadp = gr.ddp(theta)
    
    def column_int(var_in):
        var_int = mc.cp_air * var_in.sum('pfull')*5000./mc.grav * 10.**5.
        return var_int
        
    dthetadydt = column_int(gr.ddt(dthetady))
    div_term = column_int(-1.* gr.ddy(data.vcomp) * dthetady)
    horiz_adv = column_int(-1.* data.vcomp * gr.ddy(dthetady))
    tilt_term = column_int(-1.* gr.ddy(data.omega, vector=False) * dthetadp)
    vert_adv = column_int(-1.* data.omega * gr.ddy(dthetadp, vector=False))
    heating_grad = column_int(gr.ddy(heating, vector=False))
    
    #eddies = column_int(-1.*gr.ddy(gr.ddy(vcomp_theta) + gr.ddp(omega_theta))) - (div_term + horiz_adv + tilt_term + vert_adv)
    
    #eddies = dthetadydt - div_term - horiz_adv - tilt_term - vert_adv - heating_grad
    
    term_sum = div_term + horiz_adv + tilt_term + vert_adv + heating_grad
    
    # Six subplots
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f1=dthetadydt.plot.contourf(ax=ax1, x='xofyear', y='lat', extend='both', levels=np.arange(-2.,2.1,0.2), add_labels=False)
    psi.plot.contour(ax=ax1, x='xofyear', y='lat', extend='both', add_labels=False, colors='k', alpha=0.2, levels=np.arange(0,600.,100.))
    ax1.set_ylabel('Latitude')
    ax1.set_title('$d\Theta/dtdy$', fontsize=17)
    ax1.set_ylim(-60,60)
    ax1.grid(True,linestyle=':')
    ax1.set_yticks(np.arange(-60.,61.,30.))
    ax1.text(-15, 60, 'a)')
    
    #Second plot
    div_term.plot.contourf(ax=ax2, x='xofyear', y='lat', extend='both', levels=np.arange(-9.,9.1,0.5), add_labels=False)
    psi.plot.contour(ax=ax2, x='xofyear', y='lat', extend='both', add_labels=False, colors='k', alpha=0.2, levels=np.arange(0,600.,100.))
    ax2.grid(True,linestyle=':')
    ax2.set_title('$-dv/dy*d\Theta/dy$', fontsize=17)
    ax2.set_ylim(-60,60)
    ax2.text(-5, 60, 'b)')
    
    #Third plot
    horiz_adv.plot.contourf(ax=ax3, x='xofyear', y='lat', extend='both', levels=np.arange(-9.,9.1,0.5), add_labels=False)
    psi.plot.contour(ax=ax3, x='xofyear', y='lat', extend='both', add_labels=False, colors='k', alpha=0.2, levels=np.arange(0,600.,100.))
    ax3.grid(True,linestyle=':')
    ax3.set_ylim(-60,60)
    ax3.set_title('$-vd\Theta/dydy$', fontsize=17)
    ax3.text(-5, 60, 'c)')
    
    #Fourth plot
    tilt_term.plot.contourf(ax=ax4, x='xofyear', y='lat', extend='both', levels=np.arange(-9.,9.1,0.5), add_labels=False)
    psi.plot.contour(ax=ax4, x='xofyear', y='lat', extend='both', add_labels=False, colors='k', alpha=0.2, levels=np.arange(0,600.,100.))
    ax4.grid(True,linestyle=':')
    ax4.set_ylabel('Latitude')
    ax4.set_ylim(-60,60)
    ax4.set_title('$-dw/dy*d\Theta/dp$', fontsize=17)
    ax4.set_yticks(np.arange(-60.,61.,30.))
    ax4.text(-15, 60, 'd)')
    
    #Fifth plot
    vert_adv.plot.contourf(ax=ax5, x='xofyear', y='lat', extend='both', levels=np.arange(-9.,9.1,0.5), add_labels=False)
    psi.plot.contour(ax=ax5, x='xofyear', y='lat', extend='both', add_labels=False, colors='k', alpha=0.2, levels=np.arange(0,600.,100.))
    ax5.grid(True,linestyle=':')
    ax5.set_title('$-wd\Theta/dydp$', fontsize=17)
    ax5.set_ylim(-60,60)
    ax5.text(-5, 60, 'e)')
    
    #Sixth plot
    heating_grad.plot.contourf(ax=ax6, x='xofyear', y='lat', extend='both', levels=np.arange(-9.,9.1,0.5), add_labels=False)
    psi.plot.contour(ax=ax6, x='xofyear', y='lat', extend='both', add_labels=False, colors='k', alpha=0.2, levels=np.arange(0,600.,100.))
    ax6.grid(True,linestyle=':')
    ax6.set_ylim(-60,60)
    ax6.set_title('Diabatic heating gradient', fontsize=17)
    ax6.text(-5, 60, 'f)')
    
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0.1, hspace=0.25, wspace=0.12)
    #Colorbar
    #cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    #cb1=fig.colorbar(f1, ax=[ax1,ax2,ax3,ax4,ax5,ax6], use_gridspec=True,fraction=0.15, aspect=30)
    #cb1.set_label('$ms^{-1}day^{-1}$')
    
    if moist:
        figname = 'temp_grad_tend_' +run+ '_moist.pdf'
    else:
        figname = 'temp_grad_tend_' +run+ '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
temp_grad_tend('sn_1.000_evap_fluxes_heattrans')
temp_grad_tend('sn_1.000_evap_fluxes_heattrans', moist=True)
