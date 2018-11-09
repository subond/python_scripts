"""
20/06/2018
Alternative vorticity breakdown - separate planetary and relative vorticity terms
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

def vort_budg_terms(run, lonin=[-1.,361.], rot_fac=1.):
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
    
    
    # Calculate vertical component of pentad mean vorticity parts: f and dv/dx - du/dy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    
    # Evaluate relative vorticity
    vor = v_dx - u_dy 
    
    # Take divergencc of u*vor
    dvorudx = gr.ddx(vor*data.ucomp)
    dvorvdy = gr.ddy(vor*data.vcomp)
    
    div_uvor = -86400.**2 * (dvorudx + dvorvdy)
    
    dfudx = gr.ddx(data.ucomp * f)
    dfvdy = gr.ddy(data.vcomp * f)
    
    div_uf = -86400.**2 * (dfudx + dfvdy)
    

    # If run is ap_2 or full_qflux, then vorticity budget and velocities were saved separately. Load the budget up now
    if run in ['ap_2','full_qflux']:
        data = data_vort
        # Another glitch with ap_2 and full_qflux - these are only on 150 hPa.
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
    output_dict = {'div_uvor': (dim_list, div_uvor.sel(lon=lons).mean('lon')),
                 'div_uf':  (dim_list, div_uf.sel(lon=lons).mean('lon'))}
                 
    
    if run in ['ap_2','full_qflux']:
        # Single level transients from ap_2, full_qflux mistakes
        transient_hm = ((data.stretching.sel(pfull=150).values + data.horiz_md.sel(pfull=150).values) * 86400.**2.
                          - (div_uvor + div_uf).sel(pfull=150))
    else:
        # Average out longitude for other plot types 
        transient_hm = ((data.stretching + data.horiz_md)*86400.**2. - div_uvor - div_uf).sel(lon=lons).mean('lon')
    output_dict.update({'transient_hm':  (dim_list, transient_hm)})

    # Create a dataset of the terms to be plotted and return
    ds = xr.Dataset(output_dict, coords=coord_dict)

    return ds
    


def vort_eq_hm(run, lev=150., lonin=[-1.,361.], month_labels=True, rot_fac=1.,  add_psi=True):
    '''Plot vorticity budget Hovmoller. RG 3/11/2017
       Imputs: run = run_name
               lev = level to plot
               lonin = longitude range to average over
               planetary_only = only plot planetary vorticity terms
               month_labels = label x axis with month of year (no good for non earth year length)
               rot_fac = scale factor for Earth's rotation rate
               no_eddies = don't plot transient eddies too'''
    
    # Call the above function to get fields to plot. For a Hovmoller, run is by definition not steady state. Pass other inputs on.
    ds = vort_budg_terms(run, lonin=lonin, rot_fac=rot_fac)
    
    
    if add_psi:
        data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
        if lonin[1]>lonin[0]:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
        else:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        psi = mass_streamfunction(data, a=6376.0e3, dp_in=50., lons=lons)
        psi /= 1.e9
    
    
    # Determine figure size based on number of panels
    rcParams['figure.figsize'] = 15, 4
        
    rcParams['font.size'] = 14
    rcParams['text.usetex'] = True
    
    levels = np.arange(-1.5,1.6,0.25)
    
    print('starting plotting')
    
    # Set number of panels and declare which are the left column, bottom row.
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey='row')
    left_column = [ax1]
    all_plots = [ax1,ax2,ax3]
        
    plt.set_cmap('RdBu_r')
    
    # Select the pressure level, if needed
    if run in ['ap_2','full_qflux']:
        ds = ds
    else:
        ds = ds.sel(pfull=lev)
    
    # Plot! 
    
    f2 = ds.div_uvor.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax1, extend = 'both', add_labels=False)
    if add_psi:
        ax1.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=500), levels=np.arange(-500.,510.,100.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
        ax1.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=500), levels=np.arange(-500.,510.,500.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
    ax1.set_title('Relative vorticty part', fontsize=17)
    ax1.text(-15, 60, 'a)')
    
    ds.div_uf.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax2, extend = 'both', add_labels=False)
    if add_psi:
        ax2.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=500), levels=np.arange(-500.,510.,100.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
        ax2.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=500), levels=np.arange(-500.,510.,500.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
    ax2.set_title('Planetary vorticity part', fontsize=17)
    ax2.text(-7, 60, 'b)')
    
    ds.transient_hm.plot.contourf(x='xofyear', y='lat', levels=levels, ax=ax3, extend = 'both', add_labels=False)
    if add_psi:
        ax3.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=150), levels=np.arange(-500.,510.,100.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
        ax3.contour(data.xofyear, data.lat, -1.*psi.sel(pfull=150), levels=np.arange(-500.,510.,500.), add_labels=False, alpha=0.15, colors='k', linewidths=2)
    ax3.set_title('Transient eddy vorticity tendency', fontsize=17)
    ax3.text(-7, 60, 'c)')
    
    # Add grid, set y limits and ticks
    for ax in all_plots:
        ax.grid(True,linestyle=':')
        ax.set_ylim(-60,60)
        ax.set_yticks(np.arange(-60,61,30))
        if month_labels:
            mn_dic = month_dic(1)
            tickspace = list(range(13,72,18))
            labels = [mn_dic[(k+5)/6 ] for k in tickspace]
            ax.set_xticks(tickspace)
            ax.set_xticklabels(labels,rotation=25)
    
    # Label left column
    for ax in left_column:
        ax.set_ylabel('Latitude')
        
            
    plt.subplots_adjust(right=0.97, left=0.08, top=0.93, bottom=0.1, hspace=0.25, wspace=0.15)
    
    # Set plot name based on type of plot.
    if lonin == [-1.,361.]:
        lon_tag = ''
    else:
        lon_tag = '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1]))
    
        
    figname = 'alt_breakdown_' + run + lon_tag + '.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    

if __name__ == "__main__":
    
    vort_eq_hm('control_qflux', lonin=[60.,150.])
    vort_eq_hm('sn_1.000')
    

