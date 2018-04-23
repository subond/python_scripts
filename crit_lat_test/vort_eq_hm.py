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
    
    rcParams['figure.figsize'] = 10, 7.5
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/vort_schematic/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Load in vorticity budget term means
    data_vort = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_'+run+'.nc')
    data_vort = data_vort * 86400.**2. #Convert to day^-2
    
    print 'vorticity budget data loaded'
    
    #Also load climatological data so that transient eddies can be calculated 
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_uv'+run+'.nc')
    print 'climatology loaded'
    
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
    
    transient = data_vort.horiz_md.sel(pfull=lev).values + data_vort.stretching.sel(pfull=lev).values - horiz_md_mean.values - stretching_mean.values
    
    data_vort['transient'] = (('pentad','lat','lon'), transient )	
    
    horiz_md_hm = horiz_md_mean.sel(lon=lons).mean('lon')
    stretching_hm = stretching_mean.sel(lon=lons).mean('lon')
    transient_hm = data_vort.transient.sel(lon=lons).mean('lon')
    total_hm = data_vort.total.sel(pfull=lev, lon=lons).mean('lon')
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    levels = np.arange(-1.5,1.6,0.1)
    
    print 'starting plotting'

    # Four subplots
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f2=horiz_md_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax1, extend = 'both', add_labels=False, add_colorbar=False)
    ax1.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)
    ax1.set_yticks(np.arange(-60,61,30))
    ax1.grid(True,linestyle=':')
    ax1.text(-10, 60, 'a)')

    #Second plot
    stretching_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax2, extend = 'both', add_labels=False, add_colorbar=False)
    ax2.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax2.set_ylim(-60,60)
    ax2.set_yticks(np.arange(-60,61,30))
    ax2.grid(True,linestyle=':')
    ax2.text(-5, 60, 'b)')
    
    
    #Third plot
    transient_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax3, extend = 'both', add_labels=False, add_colorbar=False)
    ax3.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax3.set_ylabel('Latitude')
    ax3.set_ylim(-60,60)
    ax3.set_yticks(np.arange(-60,61,30))
    ax3.set_xticks(tickspace)
    ax3.set_xticklabels(labels,rotation=25)
    ax3.grid(True,linestyle=':')
    ax3.text(-10, 60, 'c)')
    
    #Fourth plot
    total_hm.plot.contourf(x='pentad', y='lat', levels=levels, ax=ax4, extend = 'both', add_labels=False, add_colorbar=False)
    ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax4.set_ylim(-60,60)
    ax4.set_yticks(np.arange(-60,61,30))
    ax4.set_xticks(tickspace)
    ax4.set_xticklabels(labels,rotation=25)
    ax4.grid(True,linestyle=':')
    ax4.text(-5, 60, 'd)')
    
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.2, wspace=0.1)
    #Colorbar
    cb1=f.colorbar(f2, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Vorticity tendency, day$^{-2}$')
    
    
    if lonin == [-1.,361.]:
        figname = 'vort_budg_hm_' + run + '.pdf'
    else:
        figname = 'vort_budg_hm_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
      
vort_eq_hm('ap_2')
#vort_eq_hm('full_qflux')
vort_eq_hm('full_qflux', lonin=[60.,150.])


