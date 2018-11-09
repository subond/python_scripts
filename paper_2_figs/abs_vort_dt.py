"""
Functions to look at absolute vorticity evolution in relation to precip centroid latitude (19/04/2018)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, make_sym, cell_area
from climatology import precip_centroid
from pylab import rcParams
import scipy.interpolate as spint
from hadley_cell import mass_streamfunction

def psi_mean_clim(run):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    
    area = cell_area(42, '/scratch/rg419/Isca/')
    area_xr = xr.DataArray(area, [('lat', data.lat ), ('lon', data.lon)])
    area_xr = area_xr.mean('lon')
    psi_mean = ((psi*area_xr).sum(('lat'))/area_xr.sum(('lat'))).mean('pfull')
    #psi_mean = psi_mean*-1.

    return psi_mean
    
def abs_vort_dt_plot(run, ax, rot_fac=1., lev=150., dvordt_flag=False):
    '''Plot dvordt or 1/vor * dvordt'''
    
    #Load data
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    #Calculate psi to overplot
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    psi = make_sym(psi, asym=True)
    
    # Make precip symmetric and find the precip centroid
    data['precipitation'] = make_sym(data.precipitation)
    precip_centroid(data)
    
    # Calculate vorticity
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy + f).sel(pfull=lev)*86400.

    vor = make_sym(vor, asym=True)
    
    # Take time derivative of absolute vorticity
    dvordt = gr.ddt(vor.mean('lon'))*86400.
    # Also normalise this by the value of vorticity
    dvordtvor = dvordt/vor.mean('lon')
    
    # Plot!
    plot_dir = '/scratch/rg419/plots/paper_2_figs/abs_vort_dt/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    rcParams['figure.figsize'] = 5.5, 4.3
    rcParams['font.size'] = 14

    #fig = plt.figure()
    #ax1 = fig.add_subplot(111)    
    
    if dvordt_flag:
        f1 = dvordt.plot.contourf(ax=ax, x='xofyear', y='lat', levels=np.arange(-0.06,0.07,0.01), add_colorbar=False, add_labels=False)
    else:
        f1 = dvordtvor.plot.contourf(ax=ax, x='xofyear', y='lat', levels=np.arange(-0.05,0.055,0.005),  add_colorbar=False, add_labels=False)
    #psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    #psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    #psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax, color='k', linewidth=2)
    ax.set_ylim([-60,60])
    ax.set_ylabel('Latitude')
    ax.set_xticks([12,24,36,48,60,72])
    ax.set_yticks([-60,-30,0,30,60])
    ax.set_xlabel('Pentad')
    ax.grid(True,linestyle=':')
    
    #originalSize = get(gca, 'Position')
    #set(f1, 'Position', originalSize)
    
    #plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.05)
    
    box = ax.get_position()
    #ax.set_position([box.x0*1.05, box.y0, box.width, box.height])
    
    axColor = plt.axes([box.x0 + box.width*0.92, box.y0+box.y0*0.12, 0.015, box.height])
    
    #cb1=fig.colorbar(f1, ax=ax, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.2, aspect=40)
    #cb1=plt.colorbar(f1, ax=axColor, use_gridspec=True, orientation = 'vertical',fraction=0.15, pad=0.05, aspect=30)
    cb1=fig.colorbar(f1, cax=axColor, orientation = 'vertical')
    
    #if dvordt_flag:
    #    cb1.set_label('Absolute vorticity tendency, day$^{-2}$')
        #plt.savefig(plot_dir+'abs_vordt_' + run + '.pdf', format='pdf')
    #else:
    #    cb1.set_label('Normalised absolute vorticity tendency, day$^{-1}$')
        #plt.savefig(plot_dir+'abs_vordtvor_' + run + '.pdf', format='pdf')
    #plt.close()
        

def time_interp(data, times_new):
    fun = spint.interp1d(data.xofyear.values, data.values, axis=0, fill_value='extrapolate', kind='quadratic')
    data_new = fun(times_new)
    data_new = xr.DataArray(data_new, coords=[times_new, data.lat.values, data.lon.values], dims=['xofyear','lat','lon'])
    return data_new
    
    
def lat_interp(data, lats_new):
    fun = spint.interp1d(data.lat.values, data.values, axis=-1, fill_value='extrapolate', kind='quadratic')
    data_new = fun(lats_new)
    data_new = xr.DataArray(data_new, coords=[data.xofyear.values, lats_new], dims=['xofyear','lat'])
    return data_new
    

def abs_vort_dt_at_pcent(run, rot_fac=1., lev=150.,lat_bound=45., res=0.01, interp=True):
    '''Calculate the (normalised) vorticity tendency at the precipitation centroid'''
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    # Interpolate in time if wanted
    if interp:
        times_new = np.arange(1., 72.2, 0.2)
        precip_new = time_interp(data.precipitation, times_new)
        u_new = time_interp(data.ucomp.sel(pfull=lev), times_new)
        v_new = time_interp(data.vcomp.sel(pfull=lev), times_new)
    else:
        precip_new = data.precipitation
        u_new = data.ucomp.sel(pfull=lev)
        v_new = data.vcomp.sel(pfull=lev)
        
    # Find precipitation centroid
    precip_new = make_sym(precip_new)
    data = xr.Dataset({'precipitation': precip_new}, coords=precip_new.coords)
    precip_centroid(data, lat_bound=lat_bound, res=res)
    
    # Calculate vorticity
    v_dx = gr.ddx(v_new)  # dvdx
    u_dy = gr.ddy(u_new)  # dudy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy + f)*86400.
    div = gr.ddx(u_new) + gr.ddy(v_new)
    stretching_mean = (-86400. * vor * div).mean('lon')
    vor = make_sym(vor, asym=True)
    stretching_mean = make_sym(stretching_mean, asym=True)
    
    # Take time derivative of absolute vorticity
    if interp:
        dvordt = gr.ddt(vor.mean('lon'))*86400.*5.
    else:
        dvordt = gr.ddt(vor.mean('lon'))*86400.
    # Also normalise this by the value of vorticity
    dvordtvor = dvordt/vor.mean('lon')
    #dvordtvor = stretching_mean/vor.mean('lon')
    
    # Interpolate vorticity in latitude to match precipitation centroid lats
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -lat_bound and data.lat[i] <= lat_bound]    
    lats_new = np.arange(-lat_bound, lat_bound+res, res)
    stretching_mean = lat_interp(stretching_mean.sel(lat=lats), lats_new)
    dvordt = lat_interp(dvordt.sel(lat=lats), lats_new)
    dvordtvor = lat_interp(dvordtvor.sel(lat=lats), lats_new)
    
    # Get and return dvordt and dvordtvor at the precipitation centroid, as well as the precipitation centroid itself
    st_pcent = [float(stretching_mean[i,:].sel(lat=data.p_cent.values[i]).values) for i in range(len(data.xofyear))]
    st_pcent = xr.DataArray(np.asarray(st_pcent), coords=[stretching_mean.xofyear.values], dims=['xofyear'])

    dvordt_pcent = [float(dvordt[i,:].sel(lat=data.p_cent.values[i]).values) for i in range(len(data.xofyear))]
    dvordt_pcent = xr.DataArray(np.asarray(dvordt_pcent), coords=[dvordt.xofyear.values], dims=['xofyear'])
    
    dvordtvor_pcent = [float(dvordtvor[i,:].sel(lat=data.p_cent.values[i]).values) for i in range(len(data.xofyear))]
    dvordtvor_pcent = xr.DataArray(np.asarray(dvordtvor_pcent), coords=[dvordtvor.xofyear.values], dims=['xofyear'])
    
    return dvordt_pcent, dvordtvor_pcent, data.p_cent, st_pcent


def abs_vort_at_pcent_plot(runs, rot_facs, colors, ax=None, lev=150., filename='abs_vort_dtvor_pcent', interp=True, dvordt=False, period_facs=None):
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    rcParams['figure.figsize'] = 6.5, 4
    rcParams['font.size'] = 14
    
    if ax==None:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax=ax1
    
    if period_facs==None:
        period_facs = [1.]*len(runs)
    
    for i in range(len(runs)):
        dvordt_pcent, dvordtvor_pcent, pcent, st_pcent = abs_vort_dt_at_pcent(runs[i], rot_facs[i], interp=interp)
        if interp:
            dpcentdt = gr.ddt(pcent)*86400.*5.*period_facs[i]
        else:
            dpcentdt = gr.ddt(pcent)*86400.*period_facs[i]
            
        times = [dvordt_pcent.xofyear[j] for j in range(len(dvordt_pcent.xofyear)) if dvordt_pcent.xofyear[j] >= 40*period_facs[i] and dvordt_pcent.xofyear[j] <= 55*period_facs[i]]
        #times = [dvordt_pcent.xofyear[i] for i in range(len(dvordt_pcent.xofyear)) if dvordt_pcent.xofyear[i] >= 0 and dvordt_pcent.xofyear[i] <= 10000]
        dvordtvor_pcent = dvordtvor_pcent.where(((dvordtvor_pcent < 0.02) & (dvordtvor_pcent > -0.04)))
        if dvordt:
            #ax.plot(dpcentdt.sel(xofyear=times), dvordt_pcent.sel(xofyear=times), color=colors[i], linewidth=2)
            ax.plot(pcent.sel(xofyear=times), dvordt_pcent.sel(xofyear=times), color=colors[i], linewidth=2)
            ax.set_ylabel('d$\zeta$/dt at precip. centroid')     
        else:
            ax.plot(dpcentdt.sel(xofyear=times), dvordtvor_pcent.sel(xofyear=times), color=colors[i], linewidth=2)
            #ax.plot(pcent.sel(xofyear=times), dvordtvor_pcent.sel(xofyear=times), color=colors[i], linewidth=2)
            ax.set_ylim([-0.04,0.02])
            #ax.set_ylim([-0.1,0.1])
            ax.set_ylabel('$1/(f+\zeta)$*d$\zeta$/dt at ITCZ')
    
    ax.set_xlabel('ITCZ migration rate')
    ax.set_xlim([0.,0.6])
    ax.grid(True,linestyle=':')
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    legend = ax.legend(rot_facs, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='$\Omega$/$\Omega_{E}$', fontsize=10)
    legend.get_title().set_fontsize(10)
    
    if ax==None:
        plt.subplots_adjust(left=0.18, right=0.8, top=0.95, bottom=0.15)

        plt.savefig(plot_dir+filename+'.pdf', format='pdf')
        plt.close()


if __name__ == "__main__":

    #colors=['b','g','k','r','c','m','y']
    #runs = ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']
    #rot_facs=[0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.]
    #abs_vort_at_pcent_plot(runs, rot_facs, colors)
    #abs_vort_at_pcent_plot(runs, rot_facs, colors, filename='abs_vort_dtvor_pcent_nointerp', interp=False)
    
    #colors=['b','g','k','r','c']
    #runs = ['mld_2.5', 'mld_5', 'sn_1.000', 'mld_15', 'mld_20']
    #rot_facs=[1.]*5
    #abs_vort_at_pcent_plot(runs, rot_facs, colors, filename='abs_vort_dtvor_pcent_mld', interp=False)
    
    #colors=['g','k','r','c','m']
    #runs = ['sn_0.500', 'sn_1.000', 'sn_2.000', 'sn_3.000', 'sn_4.000']
    #rot_facs=[1.]*5
    #period_facs=[0.5,1.,2.,3.,4.]
    #abs_vort_at_pcent_plot(runs, rot_facs, colors, filename='abs_vort_dtvor_pcent_sn', interp=False, period_facs=period_facs)
    
    
    
    # Set plotting directory
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)

    # Set figure parameters
    rcParams['figure.figsize'] = 5.5, 8.5
    rcParams['font.size'] = 14

    # Start figure with 4 subplots
    fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1)
    
    control = 'sine_sst_10m'#'sn_1.000'
    abs_vort_dt_plot(control, ax=ax1, dvordt_flag=True)

    psi_mean = psi_mean_clim(control)
    dvordt_pcent, dvordtvor_pcent, p_cent, st_pcent =  abs_vort_dt_at_pcent(control, interp=False)
    dpcentdt = gr.ddt(p_cent)*86400.
    dpcentdt.plot(ax=ax2, color='k', linewidth=2)
    #p_cent.plot(ax=ax2, color='k', linewidth=2)
    
    #dpsi_meandt = gr.ddt(psi_mean)*86400./5.
    #dpsi_meandt.plot(ax=ax2, color='g', linewidth=2)
    
    ax2_twin = ax2.twinx()
    #st_pcent.plot(ax=ax2_twin, color='b', linewidth=2)
    dvordt_pcent.plot(ax=ax2_twin, color='b', linewidth=2)
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax2.set_xlabel('Pentad')
    ax2.set_ylabel('ITCZ migration rate')
    ax2_twin.set_ylabel('d$\zeta$/dt at ITCZ', color='b')
    ax2.set_xlim([0,72])
    ax2.set_xticks(np.arange(0,73,12))
    #ax2.set_yticks(np.arange(-30.,30.1,10.))
    ax2_twin.set_yticks(np.arange(-0.06,0.07,0.02))
    ax2_twin.spines['right'].set_color('blue')
    ax2_twin.yaxis.label.set_color('blue')
    [t.set_color('blue') for t in ax2_twin.yaxis.get_ticklabels()]
    ax2.grid(True,linestyle=':')
    
    colors=['b','g','k','r','c','m','y']
    runs = ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']
    #runs_5 = ['rt_0.500_5', 'rt_0.750_5', 'mld_5', 'rt_1.250_5', 'rt_1.500_5', 'rt_1.750_5', 'rt_2.000_5']        
    #runs_15 = ['rt_0.500_15', 'rt_0.750_15', 'mld_15', 'rt_1.250_15', 'rt_1.500_15', 'rt_1.750_15', 'rt_2.000_15']
    
    rot_facs=[0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.]
    #rot_facs=[1.,1.,1.,1.,1.]
    #runs = ['mld_2.5', 'mld_5', 'sn_1.000', 'mld_15', 'mld_20']
    
    abs_vort_at_pcent_plot(runs, rot_facs, colors, ax=ax3, interp=False)#, dvordt=True)
    
    ax1.text(-12, 60., 'a)')
    ax2.text(-12, 0.6, 'b)')
    ax3.text(-0.12, 0.02, 'c)')
    
    plt.subplots_adjust(left=0.17, right=0.8, top=0.97, bottom=0.07, hspace=0.3)
    
    plt.savefig(plot_dir+'abs_vort_dt_' + control + '.pdf', format='pdf')
    plt.close()
    
    
    
    