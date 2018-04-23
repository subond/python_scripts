"""
Functions to look at absolute vorticity evolution in relation to precip centroid latitude (19/04/2018)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr
from climatology import precip_centroid
from pylab import rcParams
import scipy.interpolate as spint

    

def abs_vort_dt_plot(run, rot_fac=1., lev=150.):

    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    precip_centroid(data)
    
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy + f).sel(pfull=lev)*86400.
    
    # Take time derivative of absolute vorticity
    dvordt = gr.ddt(vor.mean('lon'))*86400.
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    rcParams['figure.figsize'] = 8, 4
    rcParams['font.size'] = 16

    fig = plt.figure()
    ax1 = fig.add_subplot(111)    
    
    f1 = dvordt.plot.contourf(ax=ax1, x='xofyear', y='lat', levels=np.arange(-0.06,0.07,0.01), add_colorbar=False, add_labels=False)
    data.p_cent.plot.line(color='k', linewidth=2)
    ax1.set_ylim([-60,60])
    ax1.set_ylabel('Latitude')
    ax1.set_xticks([12,24,36,48,60,72])
    ax1.set_yticks([-60,-30,0,30,60])
    ax1.set_xlabel('Pentad')
    ax1.grid(True,linestyle=':')
    
    plt.subplots_adjust(left=0.12, right=0.95, top=0.9, bottom=0.15)

    cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'vertical',fraction=0.15, pad=0.05, aspect=30)
    cb1.set_label('Absolute vorticity tendency, day$^{-2}$')

    plt.savefig(plot_dir+'abs_vort_dt_control.pdf', format='pdf')
    plt.close()
        


def abs_vort_dt_at_pcent(run, rot_fac=1., lev=150.,lat_bound=45.):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    precip_centroid(data, lat_bound=lat_bound)
    
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy + f).sel(pfull=lev)*86400.
    
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


def abs_vort_at_pcent_plot(runs, rot_facs, colors, lev=150.):
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    rcParams['figure.figsize'] = 8, 4
    rcParams['font.size'] = 16

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    for i in range(len(runs)):
        dvordt_pcent = abs_vort_dt_at_pcent(runs[i], rot_facs[i])
        dvordt_pcent.plot.line(ax=ax1, color=colors[i], linewidth=2)
    
    ax1.set_ylabel('d(vor)/dt at precip. centroid')
    ax1.set_xticks([12,24,36,48,60,72])
    ax1.set_xlabel('Pentad')
    ax1.grid(True,linestyle=':')
    
    plt.subplots_adjust(left=0.12, right=0.95, top=0.9, bottom=0.15)


    plt.savefig(plot_dir+'abs_vort_dt_pcent.pdf', format='pdf')
    plt.close()


abs_vort_dt_plot('sn_1.000')

colors=['r','m','k','c','b']
runs = ['rt_0.500', 'rt_0.750', 'sn_1.000', 
        'rt_1.500', 'rt_2.000']
rot_facs=[0.5,0.75,1.,1.5,2.]
abs_vort_at_pcent_plot(runs, rot_facs, colors)
