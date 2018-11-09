"""
3/11/2017
Contains 4 functions:
1) vort_budg_terms - evaluates pentad mean vorticity budget and transients relative to a daily snapshot budget.
2) vort_eq_hm - plots a Hovmoller plot of the vorticity budget using the above
3) vort_eq_ss - plots a latitude-pressure plot for a steady state simulation
4) vort_eq_ll - plots a latitude-longitude plot for a steady state simulation
All plots are saved in /scratch/rg419/plots/vorticity_eq_clean/
9/11/2017 - modified to have option to add psi contours to hm plots
14/02/2018 - modified to allow non steady state lat lon output
11/06/2018 - modified temporarily to read from the disca gv2 restore, and to read the qflux_ss file
"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling_updates import month_dic, gradients as gr
import sh
from pylab import rcParams
from hadley_cell import mass_streamfunction



plot_dir = '/scratch/rg419/plots/vorticity_eq_clean/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

def vort_budg_terms(run, lonin=[-1.,361.], do_ss=False, rot_fac=1., planetary_only=False, no_eddies=False, ll=False):
    '''Evaluate pentad mean vorticity budget and differences from daily snapshot budget RG 3/11/2017
       Imputs: run = run_name
               lonin = longitude range to average over
               do_ss = run is steady state
               rot_fac = scale factor for Earth's rotation rate
               planetary_only = only plot planetary vorticity terms
               no_eddies = don't plot transient eddies too
               ll = keep longitude dimension for a lat-lon plot'''
    
    
    #Load in vorticity budget term means
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/vort_eq_'+run+'.nc')
    if run in ['ap_2', 'full_qflux']:
        data_vort = data 
        data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/vort_eq_uv'+run+'.nc')
        
    print('vorticity budget data loaded')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    
    # Calculate vertical component of pentad mean absolute vorticity = f + dv/dx - du/dy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    
    # If only want to look at parts relating to planetary vorticity, substract relative part (this way turns f into a 4d array)
    if planetary_only:
        vor = v_dx - u_dy + f - v_dx + u_dy
    else:
        vor = v_dx - u_dy + f
    
    # Take gradients of vorticity
    dvordx = gr.ddx(vor)
    dvordy = gr.ddy(vor, vector=False)
    
    # Horizontal material derivative
    horiz_md_mean = -86400.**2. * (data.ucomp * dvordx + data.vcomp * dvordy)
    
    # Calculate divergence and stretching term
    u_dx = gr.ddx(data.ucomp)
    v_dy = gr.ddy(data.vcomp)
    div = u_dx + v_dy
    stretching_mean = -86400.**2. * vor * div
    
    
    if ll:
        #Keep lon dimension if want a steady state lat-lon plot
        horiz_md_av = horiz_md_mean
        stretching_av = stretching_mean
        dvordx_av = dvordx
        dvordy_av = dvordy
        vor_av = vor*86400.
        v_av = data.vcomp
        div_av = div*86400.
        u_dx_av = u_dx*86400.
        v_dy_av = v_dy*86400.
    else:
        # Take zonal mean of all over specified longitude range
        horiz_md_av = horiz_md_mean.sel(lon=lons).mean('lon')
        stretching_av = stretching_mean.sel(lon=lons).mean('lon')
        dvordx_av = dvordx.sel(lon=lons).mean('lon')
        dvordy_av = dvordy.sel(lon=lons).mean('lon')
        vor_av = vor.sel(lon=lons).mean('lon')*86400.
        v_av = data.vcomp.sel(lon=lons).mean('lon')
        div_av = div.sel(lon=lons).mean('lon')*86400.
        u_dx_av = u_dx.sel(lon=lons).mean('lon')*86400.
        v_dy_av = v_dy.sel(lon=lons).mean('lon')*86400.
    
    # If run is ap_2 or full_qflux, then vorticity budget and velocities were saved separately. Load the budget up now
    if run in ['ap_2','full_qflux']:
        data = data_vort
    
    # For a steady state lat-pressure run, only need pfull and lat
    if do_ss:
        coord_dict = {'pfull': ('pfull', data.pfull),
                        'lat': ('lat', data.lat)}
        dim_list = ['pfull', 'lat']
    
    # For a steady state lat-lon run, also need lon. Keep pfull to allow a choice of plotting level
    elif ll:
        if 'pentad' in horiz_md_av.coords:
            coord_dict = {'xofyear': ('xofyear', data.pentad),
                        'pfull': ('pfull', data.pfull),
                        'lat': ('lat', data.lat),
                        'lon': ('lon', data.lon)}
            dim_list = ['xofyear', 'pfull', 'lat', 'lon']
        else:
            coord_dict = {'pfull': ('pfull', data.pfull),
                        'lat': ('lat', data.lat),
                        'lon': ('lon', data.lon)}
            dim_list = ['pfull', 'lat', 'lon']
    
    # Another glitch with ap_2 and full_qflux - these are only on 150 hPa.
    elif run in ['ap_2','full_qflux']:
        coord_dict = {'xofyear': ('xofyear', data.pentad),
                      'lat': ('lat', data.lat)}
        dim_list = ['xofyear', 'lat']
    
    # For a Hovmoller, keep time, pressure, and lat
    else:
        coord_dict = {'xofyear': ('xofyear', data.pentad),
                        'pfull': ('pfull', data.pfull),
                          'lat': ('lat', data.lat)}
        dim_list = ['xofyear', 'pfull', 'lat']
    
    # Specify output dictionary to be written
    output_dict = {'horiz_md': (dim_list, horiz_md_av),
                 'stretching':  (dim_list, stretching_av),
                     'dvordx': (dim_list, dvordx_av),
                     'dvordy': (dim_list, dvordy_av),
                          'v': (dim_list, v_av),
                        'vor': (dim_list, vor_av),
                        'div':  (dim_list, div_av),
                        'u_dx':  (dim_list, u_dx_av),
                        'v_dy':  (dim_list, v_dy_av)}
                 
    if not (planetary_only or no_eddies):
        if run in ['ap_2','full_qflux']:
            # Single level transients from ap_2, full_qflux mistakes
            transient_s_hm = (data.stretching.sel(pfull=150).values * 86400.**2. - stretching_mean).sel(lon=lons).mean('lon')
            transient_h_hm = (data.horiz_md.sel(pfull=150).values * 86400.**2. - horiz_md_mean).sel(lon=lons).mean('lon')
        elif ll:
            # Lat lon transients
            transient_s_hm = data.stretching.values * 86400.**2. - stretching_mean
            transient_h_hm = data.horiz_md.values * 86400.**2. - horiz_md_mean
        else:
            # Average out longitude for other plot types 
            transient_s_hm = (data.stretching.values * 86400.**2. - stretching_mean).sel(lon=lons).mean('lon')
            transient_h_hm = (data.horiz_md.values * 86400.**2. - horiz_md_mean).sel(lon=lons).mean('lon')

        transient_hm = transient_s_hm + transient_h_hm        
        output_dict.update({'transient_hm':  (dim_list, transient_hm),
                         'transient_s_hm':   (dim_list, transient_s_hm),
                         'transient_h_hm':   (dim_list, transient_h_hm)})

    # Create a dataset of the terms to be plotted and return
    ds = xr.Dataset(output_dict, coords=coord_dict)

    return ds
    


def vort_eq_hm(run, lev=150., lonin=[-1.,361.], planetary_only=False, month_labels=True, rot_fac=1., no_eddies=False, add_psi=False):
    '''Plot vorticity budget Hovmoller. RG 3/11/2017
       Imputs: run = run_name
               lev = level to plot
               lonin = longitude range to average over
               planetary_only = only plot planetary vorticity terms
               month_labels = label x axis with month of year (no good for non earth year length)
               rot_fac = scale factor for Earth's rotation rate
               no_eddies = don't plot transient eddies too'''
    
    # Call the above function to get fields to plot. For a Hovmoller, run is by definition not steady state. Pass other inputs on.
    ds = vort_budg_terms(run, lonin=lonin, rot_fac=rot_fac, planetary_only=planetary_only, no_eddies=no_eddies)
    
    
    if add_psi:
        data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
        if lonin[1]>lonin[0]:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
        else:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        psi = mass_streamfunction(data, a=6376.0e3, dp_in=50., lons=lons)
        psi /= 1.e9
    
    
    # Determine figure size based on number of panels
    if planetary_only or no_eddies:
        rcParams['figure.figsize'] = 15, 6.25
    else:
        rcParams['figure.figsize'] = 15, 8
        
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    levels = np.arange(-1.5,1.6,0.25)
    
    print('starting plotting')
    
    # Set number of panels and declare which are the left column, bottom row.
    if planetary_only or no_eddies:
        f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
        left_column = [ax1,ax4]
        bottom_row = [ax4,ax5,ax6]
        all_plots = [ax1,ax2,ax3,ax4,ax5,ax6]
    else:
        f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex='col', sharey='row')
        left_column = [ax1,ax4,ax7]
        bottom_row = [ax7,ax8,ax9]
        all_plots = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
        
    plt.set_cmap('RdBu_r')
    
    # Select the pressure level, if needed
    if run in ['ap_2','full_qflux']:
        ds = ds
    else:
        ds = ds.sel(pfull=lev)
    
    # Plot! 
    
    f2 = ds.horiz_md.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax1, extend = 'both', add_labels=False)
    if add_psi:
        ax1.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=500), levels=np.arange(-500.,510.,100.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
        ax1.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=500), levels=np.arange(-500.,510.,500.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
    ax1.set_title('Horizontal advection', fontsize=17)
    ax1.text(-15, 60, 'a)')

    ds.dvordy.plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.e-10,1.05e-10,0.1e-10), ax=ax2, extend = 'both', add_labels=False)
    if planetary_only:
        ax2.set_title('$\partial f /\partial y$', fontsize=17)
    else:
        ax2.set_title('$\partial (\overline{\zeta} + f) /\partial y$', fontsize=17)
    ax2.text(-7, 60, 'b)')
    
    ds.v.plot.contourf(x='xofyear', y='lat', levels=np.arange(-12.,13,2.), ax=ax3, extend = 'both', add_labels=False)
    ax3.set_title('$\overline{v}$', fontsize=17)
    ax3.text(-7, 60, 'c)')
    
    ds.stretching.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax4, extend = 'both', add_labels=False)
    if add_psi:
        ax4.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=500), levels=np.arange(-500.,510.,100.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
        ax4.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=500), levels=np.arange(-500.,510.,500.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
    ax4.set_title('Vortex stretching', fontsize=17)
    ax4.text(-15, 60, 'd)')
    
    ds.div.plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.,1.1,0.1), ax=ax5, extend = 'both', add_labels=False)
    ax5.set_title('Divergence', fontsize=17)
    ax5.text(-7, 60, 'e)')
    
    #ds.v_dy.plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.,1.1,0.1), ax=ax5, extend = 'both', add_labels=False)
    #ax5.set_title('dv/dy', fontsize=17)
    #ax5.text(-7, 60, 'e)')
    
    ds.vor.plot.contourf(x='xofyear', y='lat', levels=np.arange(-12.,13.,2.), ax=ax6, extend = 'both', add_labels=False)
    if planetary_only:
        ax6.set_title('Planetary vorticity', fontsize=17)
    else:
        ax6.set_title('Absolute vorticity', fontsize=17)
    ax6.text(-7, 60, 'f)')
    
    # Plot transients if needed
    if not (planetary_only or no_eddies):
        ds.transient_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax7, extend = 'both', add_labels=False)
        if add_psi:
            ax7.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=150), levels=np.arange(-500.,510.,100.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
            ax7.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=150), levels=np.arange(-500.,510.,500.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
        ax7.set_title('Transient eddy vorticity tendency', fontsize=17)
        ax7.text(-15, 60, 'g)')
    
        ds.transient_h_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax8, extend = 'both', add_labels=False)
        ax8.set_title('Transient horizontal adv.', fontsize=17)
        ax8.text(-7, 60, 'h)')
    
        ds.transient_s_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax9, extend = 'both', add_labels=False)
        ax9.set_title('Transient vortex stretching', fontsize=17)
        ax9.text(-7, 60, 'i)')
    
    # Add grid, set y limits and ticks
    for ax in all_plots:
        ax.grid(True,linestyle=':')
        ax.set_ylim(-60,60)
        ax.set_yticks(np.arange(-60,61,30))
    
    # Label left column
    for ax in left_column:
        ax.set_ylabel('Latitude')
    
    # If month labels are being used, load in the month dictionary and set up labelling, then label bottom row
    if month_labels:
        mn_dic = month_dic(1)
        tickspace = list(range(13,72,18))
        labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
        for ax in bottom_row:
            ax.set_xticks(tickspace)
            ax.set_xticklabels(labels,rotation=25)
        
            
    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.1, hspace=0.25, wspace=0.15)
    
    # Set plot name based on type of plot.
    if lonin == [-1.,361.]:
        lon_tag = ''
    else:
        lon_tag = '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1]))
    
    if planetary_only:
        f_tag = '_fonly'
    else:
        f_tag = ''
        
    figname = 'vort_breakdown_' + run + lon_tag + f_tag + '.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    



def vort_eq_ss(run, lonin=[-1.,361.], planetary_only=False, rot_fac=1., no_eddies=False):
    '''Plot steady state vorticity budget in lat-pressure coordinates. RG 3/11/2017
       Inputs: run = run_name
               lonin = longitude range to average over
               planetary_only = only plot planetary vorticity terms
               rot_fac = scale factor for Earth's rotation rate
               no_eddies = don't plot transient eddies too'''
    
    # Determine figure size based on number of panels
    if planetary_only or no_eddies:
        rcParams['figure.figsize'] = 15, 6.25
    else:
        rcParams['figure.figsize'] = 15, 8

    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    levels = np.arange(-1.5,1.6,0.25)
    
    # Call vort_budg_terms to get budget. Steady state true
    ds = vort_budg_terms(run, lonin=lonin, do_ss=True, rot_fac=rot_fac, planetary_only=planetary_only, no_eddies=no_eddies)
    
    
    print('starting plotting')
    
    # Set number of panels and declare which are the left column, bottom row.
    if planetary_only or no_eddies:
        f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
        left_column = [ax1,ax4]
        bottom_row = [ax4,ax5,ax6]
        all_plots = [ax1,ax2,ax3,ax4,ax5,ax6]
    else:
        f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex='col', sharey='row')
        left_column = [ax1,ax4,ax7]
        bottom_row = [ax7,ax8,ax9]
        all_plots = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
    
    
    plt.set_cmap('RdBu_r')
    
    # Plot!
    f2=ds.horiz_md.plot.contourf(x='lat', y='pfull', levels=levels, ax=ax1, extend = 'both', add_labels=False, yincrease=False)
    ax1.set_title('Horizontal advection', fontsize=17)
    
    ds.dvordy.plot.contourf(x='lat', y='pfull', levels=np.arange(-1.e-10,1.05e-10,0.1e-10), ax=ax2, extend = 'both', add_labels=False, yincrease=False)
    if planetary_only:
        ax2.set_title('$\partial f /\partial y$', fontsize=17)
    else:
        ax2.set_title('$\partial (\overline{\zeta} + f) /\partial y$', fontsize=17)
    
    ds.v.plot.contourf(x='lat', y='pfull', levels=np.arange(-12.,13,2.), ax=ax3, extend = 'both', add_labels=False, yincrease=False)
    ax3.set_title('$\overline{v}$', fontsize=17)
    
    ds.stretching.plot.contourf(x='lat', y='pfull', levels=levels, ax=ax4, extend = 'both', add_labels=False, yincrease=False)
    ax4.set_title('Vortex stretching', fontsize=17)
    
    ds.div.plot.contourf(x='lat', y='pfull', levels=np.arange(-1.,1.1,0.1), ax=ax5, extend = 'both', add_labels=False, yincrease=False)
    ax5.set_title('Divergence', fontsize=17)
    
    ds.vor.plot.contourf(x='lat', y='pfull', levels=np.arange(-12.,13.,2.), ax=ax6, extend = 'both', add_labels=False, yincrease=False)
    if planetary_only:
        ax6.set_title('Planetary vorticity', fontsize=17)
    else:
        ax6.set_title('Absolute vorticity', fontsize=17)
            
            
    # Plot transients if needed
    if not (planetary_only or no_eddies):
        ds.transient_hm.plot.contourf(x='lat', y='pfull', levels=levels, ax=ax7, extend = 'both', add_labels=False, yincrease=False)
        ax7.set_title('Transient eddy vorticity tendency', fontsize=17)
    
        ds.transient_h_hm.plot.contourf(x='lat', y='pfull', levels=levels, ax=ax8, extend = 'both', add_labels=False, yincrease=False)
        ax8.set_title('Transient horizontal adv.', fontsize=17)
    
        ds.transient_s_hm.plot.contourf(x='lat', y='pfull', levels=levels, ax=ax9, extend = 'both', add_labels=False, yincrease=False)
        ax9.set_title('Transient vortex stretching', fontsize=17)

    # Add grid, set x limits and ticks
    for ax in all_plots:
        ax.grid(True,linestyle=':')
        ax.set_xlim(-60,60)
        ax.set_xticks(np.arange(-60,61,30))
    
    # Label left column
    for ax in left_column:
        ax.set_ylabel('Pressure, hPa')
    
    # Label bottom row
    for ax in bottom_row:
        ax.set_xlabel('Latitude')
    
    
    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.1, hspace=0.25, wspace=0.15)
    
    # Set plot name based on type of plot.
    if lonin == [-1.,361.]:
        lon_tag = ''
    else:
        lon_tag = '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1]))
    
    if planetary_only:
        f_tag = '_fonly'
    else:
        f_tag = ''
        
    figname = 'vort_breakdown_' + run + lon_tag + f_tag + '.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    




def vort_eq_ll(run, planetary_only=False, rot_fac=1., no_eddies=False, lev=150., pcent=True, add_psi=True):
    '''Plot steady state vorticity budget in lat-lon coordinates. RG 3/11/2017
       Imputs: run = run_name
               lonin = longitude range to average over
               rot_fac = scale factor for Earth's rotation rate
               no_eddies = don't plot transient eddies too'''
    
    # Determine figure size based on number of panels
    if planetary_only or no_eddies:
        rcParams['figure.figsize'] = 15, 6.25
    else:
        rcParams['figure.figsize'] = 15, 8

    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    levels = np.arange(-2.,2.1,0.25)
    
    # Call vort_budg_terms to get budget. Lat lon true
    ds = vort_budg_terms(run, ll=True, rot_fac=rot_fac, planetary_only=planetary_only, no_eddies=no_eddies)
    # Select level to plot
    ds = ds.sel(pfull=lev)
        
    
    print('starting plotting')
    
    # Set number of panels and declare which are the left column, bottom row.
    if planetary_only or no_eddies:
        f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
        left_column = [ax1,ax4]
        bottom_row = [ax4,ax5,ax6]
        all_plots = [ax1,ax2,ax3,ax4,ax5,ax6]
    else:
        f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharex='col', sharey='row')
        left_column = [ax1,ax4,ax7]
        bottom_row = [ax7,ax8,ax9]
        all_plots = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
        
    
    plt.set_cmap('RdBu_r')
    #Plot!
    f2=ds.horiz_md.plot.contourf(x='lon', y='lat', levels=levels, ax=ax1, extend = 'both', add_labels=False)
    ax1.set_title('Horizontal advection', fontsize=17)

    ds.dvordy.plot.contourf(x='lon', y='lat', levels=np.arange(-1.e-10,1.05e-10,0.1e-10), ax=ax2, extend = 'both', add_labels=False)
    if planetary_only:
        ax2.set_title('$\partial f /\partial y$', fontsize=17)
    else:
        ax2.set_title('$\partial (\overline{\zeta} + f) /\partial y$', fontsize=17)
         
    ds.v.plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13,2.), ax=ax3, extend = 'both', add_labels=False)
    ax3.set_title('$\overline{v}$', fontsize=17)

    ds.stretching.plot.contourf(x='lon', y='lat', levels=levels, ax=ax4, extend = 'both', add_labels=False)
    ax4.set_title('Vortex stretching', fontsize=17)

    ds.div.plot.contourf(x='lon', y='lat', levels=np.arange(-1.,1.1,0.1), ax=ax5, extend = 'both', add_labels=False)
    ax5.set_title('Divergence', fontsize=17)

    ds.vor.plot.contourf(x='lon', y='lat', levels=np.arange(-12.,13.,1.), ax=ax6, extend = 'both', add_labels=False)
    if planetary_only:
        ax6.set_title('Planetary vorticity', fontsize=17)
    else:
        ax6.set_title('Absolute vorticity', fontsize=17)
        
    # Plot transients if needed
    if not (planetary_only or no_eddies):
        ds.transient_hm.plot.contourf(x='lon', y='lat', levels=levels, ax=ax7, extend = 'both', add_labels=False)
        ax7.set_title('Transient eddy vorticity tendency', fontsize=17)
    
        ds.transient_h_hm.plot.contourf(x='lon', y='lat', levels=levels, ax=ax8, extend = 'both', add_labels=False)
        ax8.set_title('Transient horizontal adv.', fontsize=17)
    
        ds.transient_s_hm.plot.contourf(x='lon', y='lat', levels=levels, ax=ax9, extend = 'both', add_labels=False)
        ax9.set_title('Transient vortex stretching', fontsize=17)
    
    
    if add_psi:
        from hadley_cell import mass_streamfunction_local
        data = xr.open_dataset('/disca/restore/gv2scratch/rg419/Data_moist/climatologies/qflux_ss/' + run + '.nc')
        data['mass_sf'] = mass_streamfunction_local(data, dp_in=50.)
        for ax in all_plots:
            (data.mass_sf/1.e9).sel(pfull=500.).plot.contour(ax=ax, x='lon', y='lat', levels=np.arange(-700.,701.,100.), add_labels=False, colors='0.7', alpha=0.5)
            ax.set_xlabel('')
            ax.set_ylabel('')

    if pcent:
        from climatology import precip_centroid_ll
        data = xr.open_dataset('/disca/restore/gv2scratch/rg419/Data_moist/climatologies/qflux_ss/' + run + '.nc')
        data = precip_centroid_ll(data)
        for ax in all_plots:
            data.p_cent.plot.line(ax=ax, color='k')
            ax.set_xlabel('')
            ax.set_ylabel('')
        
        
    #Add grid
    for ax in all_plots:
        ax.grid(True,linestyle=':')

    # Label left column
    for ax in left_column:
        ax.set_ylabel('Latitude')
    
    # Label bottom row
    for ax in bottom_row:
        ax.set_xlabel('Longitude')
    
    
    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.1, hspace=0.25, wspace=0.15)
    
    # Set plot name based on type of plot.
    if planetary_only:
        f_tag = '_fonly'
    else:
        f_tag = ''
        
    figname = 'vort_breakdown_ll_' + run + f_tag + '.pdf'


    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()



    

if __name__ == "__main__":
    
    #test = vort_budg_terms('idealised_1hill', no_eddies=False, ll=True)
    
    #vort_eq_ss('ss_95.000')
    #vort_eq_ss('ss_95.000', planetary_only=True)
    #vort_eq_ss('ss_100.000')
    #vort_eq_ss('ss_100.000', planetary_only=True)
    #vort_eq_ss('ss_105.000')
    #vort_eq_ss('ss_105.000', planetary_only=True)

    #vort_eq_hm('full_qflux', lonin=[60.,150.])
    #vort_eq_hm('full_qflux', lonin=[60.,150.], planetary_only=True)
    
    vort_eq_hm('mld_20')
    #vort_eq_hm('ob_40.000')
    #vort_eq_hm('sn_1_sst_zs')


    #vort_eq_hm('half_shallow', lonin=[340,20])
    #vort_eq_hm('half_shallow', lonin=[70,110])
    #vort_eq_hm('half_shallow', lonin=[160,200])
    #vort_eq_hm('half_shallow', lonin=[250,290])
    #vort_eq_hm('rt_1.500', rot_fac=1.5)
    
    #vort_eq_hm('rt_0.500', rot_fac=0.5)
    #vort_eq_hm('rt_0.750', rot_fac=0.75)
    #vort_eq_hm('rt_1.250', rot_fac=1.25)
    #vort_eq_hm('rt_1.500', rot_fac=1.5)
    #vort_eq_hm('rt_1.750', rot_fac=1.75)
    #vort_eq_hm('rt_0.500', planetary_only=True, rot_fac=0.5)
    #vort_eq_hm('rt_2.000', rot_fac=2.)
    #vort_eq_hm('rt_2.000', planetary_only=True, rot_fac=2.)

    #vort_eq_hm('sine_sst_10m')
    #vort_eq_hm('sn_0.500', planetary_only=True, month_labels=False)
    #vort_eq_hm('dry_zs')
    #vort_eq_hm('sn_2.000', planetary_only=True, month_labels=False)
    

