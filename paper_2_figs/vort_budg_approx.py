"""
Approximate vorticity budget while GV2 out of order
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

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    
    data['precipitation'] = make_sym(data.precipitation)
    
    precip_centroid(data)
    
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy + f).sel(pfull=lev)*86400.
    
    # Take gradients of vorticity
    dvordx = gr.ddx(vor)
    dvordy = gr.ddy(vor, vector=False)
    
    # Horizontal material derivative
    horiz_md_mean = -86400. * (data.ucomp.sel(pfull=lev) * dvordx + data.vcomp.sel(pfull=lev) * dvordy)
    
    # Calculate divergence and stretching term
    div = gr.ddx(data.ucomp.sel(pfull=lev)) + gr.ddy(data.vcomp.sel(pfull=lev))
    stretching_mean = -86400. * vor * div
    
    #vor = make_sym(vor, asym=True)
    #psi = make_sym(psi, asym=True)
    
    # Take time derivative of absolute vorticity
    dvordt = gr.ddt(vor.mean('lon'))*86400.
    stretching_mean = stretching_mean.mean('lon')
    horiz_md_mean = horiz_md_mean.mean('lon')
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/abs_vort_dt/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    rcParams['figure.figsize'] = 5, 12
    rcParams['font.size'] = 14


    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1)       
    
    f1 = dvordt.plot.contourf(ax=ax1, x='xofyear', y='lat', levels=np.arange(-0.06,0.07,0.01), add_colorbar=False, add_labels=False)
    psi.sel(pfull=500).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax1, color='k', linewidth=2)
    ax1.set_ylim([-60,60])
    ax1.set_ylabel('Latitude')
    ax1.set_xticks([12,24,36,48,60,72])
    ax1.set_yticks([-60,-30,0,30,60])
    ax1.set_xlabel('Pentad')
    ax1.grid(True,linestyle=':')
    
    f1 = stretching_mean.plot.contourf(ax=ax2, x='xofyear', y='lat', levels=np.arange(-1.5,1.6,0.25), add_colorbar=False, add_labels=False)
    psi.sel(pfull=500).plot.contour(ax=ax2, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax2, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax2, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax2, color='k', linewidth=2)
    ax2.set_ylim([-60,60])
    ax2.set_ylabel('Latitude')
    ax2.set_xticks([12,24,36,48,60,72])
    ax2.set_yticks([-60,-30,0,30,60])
    ax2.set_xlabel('Pentad')
    ax2.grid(True,linestyle=':')
    
    f1 = horiz_md_mean.plot.contourf(ax=ax3, x='xofyear', y='lat', levels=np.arange(-1.5,1.6,0.25), add_colorbar=False, add_labels=False)
    psi.sel(pfull=500).plot.contour(ax=ax3, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax3, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax3, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax3, color='k', linewidth=2)
    ax3.set_ylim([-60,60])
    ax3.set_ylabel('Latitude')
    ax3.set_xticks([12,24,36,48,60,72])
    ax3.set_yticks([-60,-30,0,30,60])
    ax3.set_xlabel('Pentad')
    ax3.grid(True,linestyle=':')
    
    f1 = (stretching_mean+horiz_md_mean).plot.contourf(ax=ax4, x='xofyear', y='lat', levels=np.arange(-1.5,1.6,0.25), add_colorbar=False, add_labels=False)
    psi.sel(pfull=500).plot.contour(ax=ax4, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax4, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax4, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax4, color='k', linewidth=2)
    ax4.set_ylim([-60,60])
    ax4.set_ylabel('Latitude')
    ax4.set_xticks([12,24,36,48,60,72])
    ax4.set_yticks([-60,-30,0,30,60])
    ax4.set_xlabel('Pentad')
    ax4.grid(True,linestyle=':')
    
    
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.05)

    cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.2, aspect=40)
    cb1.set_label('Absolute vorticity tendency, day$^{-2}$')

    plt.savefig(plot_dir+'vort_terms_' + run + '.pdf', format='pdf')
    plt.close()
        

runs = ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']
rotfacs = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.]
for i in range(7): 
    abs_vort_dt_plot(runs[i], rot_fac=rotfacs[i])
