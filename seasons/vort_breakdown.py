"""
Load in the vorticity budget terms and produce a lat-time plot at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams

def vort_eq_hm(run, lev=150, lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 15, 6.25
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/seasons/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Load in vorticity budget term means
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_'+run+'.nc')
    print 'vorticity budget data loaded'
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    
    abs_vort = vor.sel(lon=lons).mean('lon')*86400.
    
    dvordx = gr.ddx(vor)
    dvordy = gr.ddy(vor, vector=False)
    
    horiz_md_mean = -86400.**2. * (data.ucomp * dvordx + data.vcomp * dvordy)
    
    div = gr.ddx(data.ucomp) + gr.ddy(data.vcomp)
    stretching_mean = -86400.**2. * vor * div
    
    horiz_md_hm = horiz_md_mean.sel(lon=lons).mean('lon')
    stretching_hm = stretching_mean.sel(lon=lons).mean('lon')
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    levels = np.arange(-1.5,1.6,0.25)
    
    print 'starting plotting'

    # Four subplots
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f2=horiz_md_hm.sel(pfull=lev).plot.contourf(x='pentad', y='lat', levels=levels, ax=ax1, extend = 'both', add_labels=False)
    #ax1.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)
    ax1.set_yticks(np.arange(-60,61,30))
    ax1.grid(True,linestyle=':')
    ax1.text(-15, 60, 'a)')

    #Second plot
    dvordy.sel(lon=lons).mean('lon').sel(pfull=lev).plot.contourf(x='pentad', y='lat', levels=np.arange(-1.e-10,1.05e-10,0.1e-10), ax=ax2, extend = 'both', add_labels=False)
    #ax2.set_ylim(-60,60)
    ax2.set_yticks(np.arange(-60,61,30))
    ax2.grid(True,linestyle=':')
    ax2.text(-7, 60, 'b)')
    
    
    #Third plot
    data.vcomp.sel(lon=lons).mean('lon').sel(pfull=lev).plot.contourf(x='pentad', y='lat', levels=np.arange(-12.,13,2.), ax=ax3, extend = 'both', add_labels=False)
    ax3.set_ylim(-60,60)
    ax3.set_yticks(np.arange(-60,61,30))
    ax3.grid(True,linestyle=':')
    ax3.text(-7, 60, 'c)')
    
    #Fourth plot
    stretching_hm.sel(pfull=lev).plot.contourf(x='pentad', y='lat', levels=levels, ax=ax4, extend = 'both', add_labels=False)
    #ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax4.set_ylabel('Latitude')
    ax4.set_ylim(-60,60)
    ax4.set_yticks(np.arange(-60,61,30))
    #ax4.set_xticks(tickspace)
    #ax4.set_xticklabels(labels,rotation=25)
    ax4.grid(True,linestyle=':')
    ax4.text(-15, 60, 'd)')
    
    #Fith plot
    (div.sel(lon=lons).mean('lon').sel(pfull=lev)*86400.).plot.contourf(x='pentad', y='lat', levels=np.arange(-0.6,0.7,0.1), ax=ax5, extend = 'both', add_labels=False)
    ax5.set_ylim(-60,60)
    ax5.set_yticks(np.arange(-60,61,30))
   # ax5.set_xticks(tickspace)
    #ax5.set_xticklabels(labels,rotation=25)
    ax5.grid(True,linestyle=':')
    ax5.text(-7, 60, 'e)')
    
    #Sixth plot
    (vor.sel(lon=lons).mean('lon').sel(pfull=lev)*86400.).plot.contourf(x='pentad', y='lat', levels=np.arange(-12.,13.,2.), ax=ax6, extend = 'both', add_labels=False)
    ax6.set_ylim(-60,60)
    ax6.set_yticks(np.arange(-60,61,30))
    #ax6.set_xticks(tickspace)
    #ax6.set_xticklabels(labels,rotation=25)
    ax6.grid(True,linestyle=':')
    ax6.text(-7, 60, 'f)')
    
    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.1, hspace=0.25, wspace=0.15)
    
    
    if lonin == [-1.,361.]:
        figname = 'vort_breakdown_' + run + '.pdf'
    else:
        figname = 'vort_breakdown_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
      
vort_eq_hm('sn_1.000')
vort_eq_hm('sn_2.000')
vort_eq_hm('sn_0.500')


