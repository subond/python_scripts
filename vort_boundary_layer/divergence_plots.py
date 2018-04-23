"""
Plot the lower and upper level divergence for the aquaplanet time evolving runs. Is the upper level strong on eq, weak off, lower level weak on eq, strong off true for deeper MLDs?

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams

def div_hm(run, lev=900, lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 15, 6.25
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/vort_boundary_layer/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Also load climatological data so that transient eddies can be calculated 
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    print 'climatology loaded'    
    
    div = gr.ddx(data.ucomp.sel(pfull=lev)) + gr.ddy(data.vcomp.sel(pfull=lev))
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]

    
    print 'starting plotting'

    # Four subplots
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f2=horiz_md_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax1, extend = 'both', add_labels=False)
    #ax1.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)
    ax1.set_title('Horizontal advection', fontsize=17)
    ax1.set_yticks(np.arange(-60,61,30))
    ax1.grid(True,linestyle=':')
    ax1.text(-15, 60, 'a)')

    #Second plot
    dvordy.sel(lon=lons).mean('lon').plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.e-10,1.05e-10,0.1e-10), ax=ax2, extend = 'both', add_labels=False)
    #ax2.set_ylim(-60,60)
    ax2.set_yticks(np.arange(-60,61,30))
    ax2.set_title('$\partial (\overline{\zeta} + f) /\partial y$', fontsize=17)
    ax2.grid(True,linestyle=':')
    ax2.text(-7, 60, 'b)')
    
    
    #Third plot
    data.vcomp.sel(pfull=lev).sel(lon=lons).mean('lon').plot.contourf(x='xofyear', y='lat', levels=np.arange(-12.,13,2.), ax=ax3, extend = 'both', add_labels=False)
    ax3.set_ylim(-60,60)
    ax3.set_yticks(np.arange(-60,61,30))
    ax3.set_title('$\overline{v}$', fontsize=17)
    ax3.grid(True,linestyle=':')
    ax3.text(-7, 60, 'c)')
    
    #Fourth plot
    stretching_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax4, extend = 'both', add_labels=False)
    #ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax4.set_ylabel('Latitude')
    ax4.set_ylim(-60,60)
    ax4.set_yticks(np.arange(-60,61,30))
    ax4.set_xticks(tickspace)
    ax4.set_xticklabels(labels,rotation=25)
    ax4.set_title('Vortex stretching', fontsize=17)
    ax4.grid(True,linestyle=':')
    ax4.text(-15, 60, 'd)')
    
    #Fith plot
    (div.sel(lon=lons).mean('lon')*86400.).plot.contourf(x='xofyear', y='lat', levels=np.arange(-0.6,0.7,0.1), ax=ax5, extend = 'both', add_labels=False)
    ax5.set_ylim(-60,60)
    ax5.set_yticks(np.arange(-60,61,30))
    ax5.set_xticks(tickspace)
    ax5.set_xticklabels(labels,rotation=25)
    ax5.grid(True,linestyle=':')
    ax5.set_title('Divergence', fontsize=17)
    ax5.text(-7, 60, 'e)')
    
    #Sixth plot
    (vor.sel(lon=lons).mean('lon')*86400.).plot.contourf(x='xofyear', y='lat', levels=np.arange(-12.,13.,2.), ax=ax6, extend = 'both', add_labels=False)
    ax6.set_ylim(-60,60)
    ax6.set_yticks(np.arange(-60,61,30))
    ax6.set_xticks(tickspace)
    ax6.set_xticklabels(labels,rotation=25)
    ax6.grid(True,linestyle=':')
    ax6.set_title('Absolute vorticity', fontsize=17)
    ax6.text(-7, 60, 'f)')
    
    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.1, hspace=0.25, wspace=0.15)
    
    
    if lonin == [-1.,361.]:
        figname = 'vort_breakdown_b_' + run + '.pdf'
    else:
        figname = 'vort_breakdown_b_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
      
vort_eq_hm('ap_2')
vort_eq_hm('full_qflux', lonin=[60.,150.])


