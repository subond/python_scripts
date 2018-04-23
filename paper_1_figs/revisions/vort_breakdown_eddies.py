"""
Load in the vorticity budget terms and produce a lat-time plot at 150 hPa. Include eddies

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams
from hadley_cell import mass_streamfunction

def vort_eq_hm(run, lev=150, lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 15, 8
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/paper_1_figs/revisions/'
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
    vor = v_dx - u_dy + f #- v_dx + u_dy
    
    abs_vort = vor.sel(lon=lons).mean('lon')*86400.
    
    dvordx = gr.ddx(vor)
    dvordy = gr.ddy(vor, vector=False)
    
    horiz_md_mean = -86400.**2. * (data.ucomp * dvordx + data.vcomp * dvordy)
    
    div = gr.ddx(data.ucomp) + gr.ddy(data.vcomp)
    stretching_mean = -86400.**2. * vor * div
    
    horiz_md_hm = horiz_md_mean.sel(lon=lons).mean('lon')
    stretching_hm = stretching_mean.sel(lon=lons).mean('lon')
    
    # find where min occurs
    #blob_loc = stretching_hm[30:48,:].argmin()
    #test = np.unravel_index(blob_loc, np.shape(stretching_hm[30:48,:].values))
    #print stretching_hm[30+test[0],test[1]]
    
    transient_s_hm = (data_vort.stretching.sel(pfull=lev).values - stretching_mean).sel(lon=lons).mean('lon')
    transient_h_hm = (data_vort.horiz_md.sel(pfull=lev).values - horiz_md_mean).sel(lon=lons).mean('lon')
    
    transient_hm = transient_s_hm + transient_h_hm
    
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    levels = np.arange(-1.5,1.6,0.25)
    
    print 'starting plotting'
    
    data_psi = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    psi = mass_streamfunction(data_psi, a=6376.0e3, dp_in=50., lons=lons)
    psi /= 1.e9

    # Four subplots
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f2=horiz_md_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax1, extend = 'both', add_labels=False)
    ax1.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=150), levels=np.arange(-500.,510.,100.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
    ax1.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=150), levels=np.arange(-500.,510.,500.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
    #psi.sel(pfull=150).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(-500.,510.,100.), add_labels=False, alpha=0.25, colors='k', linewidths=2)
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
    data.vcomp.sel(lon=lons).mean('lon').plot.contourf(x='xofyear', y='lat', levels=np.arange(-12.,13,2.), ax=ax3, extend = 'both', add_labels=False)
    ax3.set_ylim(-60,60)
    ax3.set_yticks(np.arange(-60,61,30))
    ax3.set_title('$\overline{v}$', fontsize=17)
    ax3.grid(True,linestyle=':')
    ax3.text(-7, 60, 'c)')
    
    #Fourth plot
    stretching_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax4, extend = 'both', add_labels=False)
    ax4.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=150), levels=np.arange(-500.,510.,100.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
    ax4.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=150), levels=np.arange(-500.,510.,500.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
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
    
    #Seventh plot
    transient_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax7, extend = 'both', add_labels=False)
    ax7.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=150), levels=np.arange(-500.,510.,100.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
    ax7.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=150), levels=np.arange(-500.,510.,500.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
    #ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax7.set_ylabel('Latitude')
    ax7.set_ylim(-60,60)
    ax7.set_yticks(np.arange(-60,61,30))
    ax7.set_xticks(tickspace)
    ax7.set_xticklabels(labels,rotation=25)
    ax7.set_title('Transient eddy vorticity tendency', fontsize=17)
    ax7.grid(True,linestyle=':')
    ax7.text(-15, 60, 'g)')
    
    #Seventh plot
    transient_h_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax8, extend = 'both', add_labels=False)
    #ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax8.set_ylim(-60,60)
    ax8.set_yticks(np.arange(-60,61,30))
    ax8.set_xticks(tickspace)
    ax8.set_xticklabels(labels,rotation=25)
    ax8.set_title('Transient horizontal adv.', fontsize=17)
    ax8.grid(True,linestyle=':')
    ax8.text(-7, 60, 'h)')
    
    #Seventh plot
    transient_s_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax9, extend = 'both', add_labels=False)
    #ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax9.set_ylim(-60,60)
    ax9.set_yticks(np.arange(-60,61,30))
    ax9.set_xticks(tickspace)
    ax9.set_xticklabels(labels,rotation=25)
    ax9.set_title('Transient vortex stretching', fontsize=17)
    ax9.grid(True,linestyle=':')
    ax9.text(-7, 60, 'i)')
    
    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.1, hspace=0.25, wspace=0.15)
    
    
    if lonin == [-1.,361.]:
        figname = 'vort_breakdown_' + run + '.pdf'
    else:
        figname = 'vort_breakdown_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
      
vort_eq_hm('ap_2')
vort_eq_hm('full_qflux', lonin=[60.,150.])


