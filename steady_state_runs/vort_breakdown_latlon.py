"""
Load in the vorticity budget terms and produce a lat-time plot at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr, model_constants as mc
from pylab import rcParams

def vort_eq_ss(run):
    
    rcParams['figure.figsize'] = 15, 6.25
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/steady_state_runs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Load in vorticity budget term means
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_'+run+'.nc')
    print 'vorticity budget data loaded'
        
    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    f = 2 * mc.omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    
    abs_vort = vor*86400.
    
    dvordx = gr.ddx(vor)
    dvordy = gr.ddy(vor, vector=False)
    
    horiz_md_mean = -86400.**2. * (data.ucomp * dvordx + data.vcomp * dvordy)
    
    div = gr.ddx(data.ucomp) + gr.ddy(data.vcomp)
    stretching_mean = -86400.**2. * vor * div
    
    horiz_md = horiz_md_mean.sel(pfull=150.)
    stretching = stretching_mean.sel(pfull=150.)
        
    print 'starting plotting'
    
    levels = np.arange(-2.,2.1,0.25)
    
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f2=horiz_md.plot.contourf(x='lon', y='lat', levels=levels, ax=ax1, extend = 'both', add_labels=False)
    #ax1.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax1.set_ylabel('Latitude')
    ax1.grid(True,linestyle=':')
    #ax1.text(-15, 60, 'a)')

    #Second plot
    dvordy.sel(pfull=150.).plot.contourf(x='lon', y='lat', levels=np.arange(-1.e-10,1.05e-10,0.1e-10), ax=ax2, extend = 'both', add_labels=False)
    ax2.grid(True,linestyle=':')
    #ax2.text(-7, 60, 'b)')
    
    
    #Third plot
    data.vcomp.sel(pfull=150.).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13,2.), ax=ax3, extend = 'both', add_labels=False)
    ax3.grid(True,linestyle=':')
    #ax3.text(-7, 60, 'c)')
    
    #Fourth plot
    stretching.plot.contourf(x='lon', y='lat', levels=levels, ax=ax4, extend = 'both', add_labels=False)
    #ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax4.set_ylabel('Latitude')
    #ax4.set_xticks(tickspace)
    #ax4.set_xticklabels(labels,rotation=25)
    ax4.grid(True,linestyle=':')
    #ax4.text(-15, 60, 'd)')
    
    #Fifth plot
    (div.sel(pfull=150.)*86400.).plot.contourf(x='lon', y='lat', levels=np.arange(-1.,1.1,0.1), ax=ax5, extend = 'both', add_labels=False)
   # ax5.set_xticks(tickspace)
    #ax5.set_xticklabels(labels,rotation=25)
    ax5.grid(True,linestyle=':')
    #ax5.text(-7, 60, 'e)')
    
    #Sixth plot
    (vor.sel(pfull=150.)*86400.).plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,2.), ax=ax6, extend = 'both', add_labels=False)
    #ax6.set_xticks(tickspace)
    #ax6.set_xticklabels(labels,rotation=25)
    ax6.grid(True,linestyle=':')
    #ax6.text(-7, 60, 'f)')
    
    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.1, hspace=0.25, wspace=0.15)
    
    
    figname = 'vort_breakdown_ll_' + run + '.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()

vort_eq_ss('qflux_0_100')
vort_eq_ss('qflux_5_100')
vort_eq_ss('qflux_10_100')
vort_eq_ss('qflux_15_100')
vort_eq_ss('qflux_20_100')
vort_eq_ss('qflux_25_100')

vort_eq_ss('qflux_0_200')
vort_eq_ss('qflux_5_200')
vort_eq_ss('qflux_10_200')
vort_eq_ss('qflux_15_200')
vort_eq_ss('qflux_20_200')
vort_eq_ss('qflux_25_200')

vort_eq_ss('qflux_0_300')
vort_eq_ss('qflux_5_300')
vort_eq_ss('qflux_10_300')
vort_eq_ss('qflux_15_300')
vort_eq_ss('qflux_20_300')
vort_eq_ss('qflux_25_300')




