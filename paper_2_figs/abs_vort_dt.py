"""
Functions to look at absolute vorticity evolution in relation to precip centroid latitude (19/04/2018)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, make_sym
from climatology import precip_centroid
from pylab import rcParams
import scipy.interpolate as spint
from hadley_cell import mass_streamfunction


def abs_vort_dt_plot(run, rot_fac=1., lev=150.):

    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    
    data['precipitation'] = make_sym(data.precipitation)
    
    precip_centroid(data)
    
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy + f).sel(pfull=lev)*86400.
    
    vor = make_sym(vor, asym=True)
    psi = make_sym(psi, asym=True)
    
    # Take time derivative of absolute vorticity
    dvordt = gr.ddt(vor.mean('lon'))*86400.
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/abs_vort_dt/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    rcParams['figure.figsize'] = 5.5, 4.3
    rcParams['font.size'] = 16

    fig = plt.figure()
    ax1 = fig.add_subplot(111)    
    
    f1 = dvordt.plot.contourf(ax=ax1, x='xofyear', y='lat', levels=np.arange(-0.06,0.07,0.01), add_colorbar=False, add_labels=False)
    psi.sel(pfull=500).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(color='k', linewidth=2)
    ax1.set_ylim([-60,60])
    ax1.set_ylabel('Latitude')
    ax1.set_xticks([12,24,36,48,60,72])
    ax1.set_yticks([-60,-30,0,30,60])
    ax1.set_xlabel('Pentad')
    ax1.grid(True,linestyle=':')
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.05)

    cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.2, aspect=40)
    cb1.set_label('Absolute vorticity tendency, day$^{-2}$')

    plt.savefig(plot_dir+'abs_vort_dt_' + run + '.pdf', format='pdf')
    plt.close()
        


def abs_vort_dt_at_pcent(run, rot_fac=1., lev=150.,lat_bound=45.):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['precipitation'] = make_sym(data.precipitation)
    
    precip_centroid(data, lat_bound=lat_bound)
    
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy + f).sel(pfull=lev)*86400.
    
    vor = make_sym(vor, asym=True)
    
    # Take time derivative of absolute vorticity
    dvordt = gr.ddt(vor.mean('lon'))*86400.
    
    # Select latitudes over which to evaluate precip centroid
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -lat_bound and data.lat[i] <= lat_bound]
    dvordt_lats = dvordt.sel(lat=lats).values
    #f_lats = f.sel(lat=lats).values
    
    # Interpolate vorticity tendency and f in latitude    
    lats_new = np.arange(-lat_bound, lat_bound+0.1, 0.1)
    
    fun = spint.interp1d(lats, dvordt_lats, axis=-1, fill_value='extrapolate', kind='quadratic')
    dvordt_new = fun(lats_new)
    dvordt = xr.DataArray(dvordt_new, coords=[data.xofyear.values, lats_new], dims=['xofyear','lat'])
    
    #fun = spint.interp1d(lats, f_lats, axis=-1, fill_value='extrapolate', kind='quadratic')
    #f_new = fun(lats_new)
    #f = xr.DataArray(f_new, coords=[lats_new], dims=['lat'])
    
    dvordt_pcent = [float(dvordt[i,:].sel(lat=data.p_cent.values[i]).values) for i in range(len(data.xofyear))]
    #f_pcent = [float(f.sel(lat=data.p_cent.values[i]).values*86400.) for i in range(len(data.xofyear))]
    dvordt_pcent = xr.DataArray(np.asarray(dvordt_pcent), coords=[data.xofyear.values], dims=['xofyear'])
    
    return dvordt_pcent


def abs_vort_at_pcent_plot(runs, rot_facs, colors, lev=150., filename='abs_vort_dt_pcent'):
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    rcParams['figure.figsize'] = 6.5, 4
    rcParams['font.size'] = 16

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    for i in range(len(runs)):
        dvordt_pcent = abs_vort_dt_at_pcent(runs[i], rot_facs[i])
        dvordt_pcent.plot.line(ax=ax1, color=colors[i], linewidth=2)
    
    ax1.set_ylabel('d$\zeta$/dt at precip. centroid')
    ax1.set_xticks([12,24,36,48,60,72])
    ax1.set_xlabel('Pentad')
    ax1.grid(True,linestyle=':')
    
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    ax1.legend(rot_facs, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title='$\Omega$/$\Omega_{E}$', fontsize=14)
    
    plt.subplots_adjust(left=0.18, right=0.8, top=0.95, bottom=0.15)


    plt.savefig(plot_dir+filename+'.pdf', format='pdf')
    plt.close()


#for run in ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000', 
#            'rt_0.750_5','mld_5', 'rt_1.250_5', 'rt_0.750_15','mld_15', 'rt_1.250_15']:
#    abs_vort_dt_plot(run)

colors=['b','g','k','r','c','m','y']
runs = ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']
rot_facs=[0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.]
abs_vort_at_pcent_plot(runs, rot_facs, colors)

colors=['g','k','r']
runs = ['rt_0.750_5','mld_5', 'rt_1.250_5']
rot_facs=[0.75, 1., 1.25]    
abs_vort_at_pcent_plot(runs, rot_facs, colors, filename='abs_vort_dt_pcent_5')

runs = ['rt_0.750_15', 'mld_15', 'rt_1.250_15']
rot_facs=[0.75, 1., 1.25]
abs_vort_at_pcent_plot(runs, rot_facs, colors, filename='abs_vort_dt_pcent_15')
