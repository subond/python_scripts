'''01/11/2018 Plot up terms contributing to heat budget in flux form at a given level for an aquaplanet'''

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling_updates import gradients as gr, model_constants as mc
import sh
from pylab import rcParams
from climatology import precip_centroid
from hadley_cell import mass_streamfunction
    
    
def heat_budg_hm_fluxform(run, lev=850, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.], plot_precip=True, do_theta=False):
    
    rcParams['figure.figsize'] = 12, 8
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/heat_budg/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/'+run+'.nc')
    
    uT_dx = -86400. * gr.ddx(data.ucomp_temp.sel(pfull=lev)).mean('lon')
    vT_dy = -86400. * gr.ddy(data.vcomp_temp.sel(pfull=lev)).mean('lon')
    wT_dp = -86400. * gr.ddp(data.omega_temp).sel(pfull=lev).mean('lon')
    
    uT_eddy = data.ucomp_temp - data.ucomp * data.temp
    vT_eddy = data.vcomp_temp - data.vcomp * data.temp
    wT_eddy = data.omega_temp - data.omega * data.temp

    #uT_dx = -86400. * (data.ucomp * gr.ddx(data.temp)).sel(pfull=lev).mean('lon')
    #vT_dy = -86400. * (data.vcomp * gr.ddy(data.temp, vector=False)).sel(pfull=lev).mean('lon')
    #wT_dp = -86400. * (data.omega * gr.ddp(data.temp)).sel(pfull=lev).mean('lon')
    
    uT_dx_eddy = -86400. * gr.ddx(uT_eddy.sel(pfull=lev)).mean('lon')
    vT_dy_eddy = -86400. * gr.ddy(vT_eddy.sel(pfull=lev)).mean('lon')
    wT_dp_eddy = -86400. * gr.ddp(wT_eddy).sel(pfull=lev).mean('lon')
    
    tdt_rad = data.tdt_rad.sel(pfull=lev).mean('lon')*86400.
    tdt_conv = data.dt_tg_convection.sel(pfull=lev).mean('lon')*86400.
    tdt_cond = data.dt_tg_condensation.sel(pfull=lev).mean('lon')*86400.
    tdt_diff = data.dt_tg_diffusion.sel(pfull=lev).mean('lon')*86400.
    
    diabatic = tdt_rad + tdt_conv + tdt_cond + tdt_diff
    eddies = uT_dx_eddy + vT_dy_eddy + wT_dp_eddy
    
    rho = data.pfull*100./data.temp/mc.rdgas
    asc_cool = (data.omega /(mc.cp_air * rho)).sel(pfull=lev).mean('lon') *86400.
    heat_sum = uT_dx + vT_dy + wT_dp + diabatic + asc_cool #+ eddies
    Tdt = gr.ddt(data.temp.sel(pfull=lev)).mean('lon') *86400.
    
    horiz_adv = uT_dx + vT_dy
    vertical_term = wT_dp + asc_cool
    
    levels = np.arange(-100,101.1,10.)
    #levels = np.arange(-20,21.,2.)
    #levels = np.arange(-3,3.1,0.2)
    #levels = np.arange(-1,1.1,0.1)
    
    # Nine subplots
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    
    plot_vars = [horiz_adv, vertical_term, diabatic, heat_sum, Tdt, eddies]
                 
    axes = [ax1,ax2,ax3,ax4,ax5,ax6]
    labels = ['a)','b)','c)','d)','e)','f)']
    for i in range(6):
        f1=plot_vars[i].plot.contourf(ax=axes[i], x='xofyear', y='lat', extend='both', add_labels=False, cmap='RdBu_r', levels = levels, add_colorbar=False)
        axes[i].set_ylim(-60,60)
        axes[i].grid(True,linestyle=':')
        axes[i].set_yticks(np.arange(-60.,61.,30.))
    
    for i in [0,3]:
        axes[i].set_ylabel('Latitude')
        axes[i].text(-15, 60, labels[i])
    
    for i in [1,2,4,5]:
        axes[i].text(-5, 60, labels[i])
    
    #for i in [3,4,5]:
    #    axes[i].set_xticklabels(ticklabels,rotation=25)
        
    # set titles    
    ax1.set_title('Horizontal advection', fontsize=17)
    ax2.set_title('Vertical advection and expansion', fontsize=17)
    ax3.set_title('Diabatic', fontsize=17)
    ax4.set_title('Residual', fontsize=17)
    ax5.set_title('Temperature tendency', fontsize=17)
    ax6.set_title('Eddy convergence', fontsize=17)
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.25, wspace=0.12)
    #Colorbar
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    
    figname = 'heat_budg_fluxform_' +run+ '.pdf'
        
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
   
if __name__ == "__main__":
    
    heat_budg_hm_fluxform('sn_1.000_evap_fluxes_heattrans',lev=850)
    heat_budg_hm_fluxform('ap_2',lev=850)

