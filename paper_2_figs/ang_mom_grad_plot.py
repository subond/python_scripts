"""
Plot f - du/dy - w/v du/dp : actual criteria for AM conserving flow (8/11/2018)
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


def ang_mom_grad_plot(run, rot_fac=1., lev=150.):
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
    data['ucomp'] = make_sym(data.ucomp)
    data['vcomp'] = make_sym(data.vcomp, asym=True)
    data['omega'] = make_sym(data.omega)
    
    # Calculate vorticity
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    u_dp = gr.ddp(data.ucomp)  # dudp
    
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy + f).sel(pfull=lev).mean('lon')*86400.
    
    vertical_term = (-1.*data.omega/data.vcomp * u_dp).sel(pfull=lev).mean('lon')*86400.
        
    ang_mom_grad = vor + vertical_term
    
    # Plot!
    plot_dir = '/scratch/rg419/plots/paper_2_figs/ang_mom_grad/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    rcParams['figure.figsize'] = 5.5, 4.3
    rcParams['font.size'] = 14

    fig = plt.figure()
    ax1 = fig.add_subplot(111)    
    
    f1 = ang_mom_grad.plot.contourf(ax=ax1, x='xofyear', y='lat', levels=np.arange(-6.,6.1,1.), add_colorbar=False, add_labels=False, extend='both')
    psi.sel(pfull=lev).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=lev).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=lev).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax1, color='k', linewidth=2)
    ax1.set_ylim([-60,60])
    ax1.set_ylabel('Latitude')
    ax1.set_xticks([12,24,36,48,60,72])
    ax1.set_yticks([-60,-30,0,30,60])
    ax1.set_xlabel('Pentad')
    ax1.grid(True,linestyle=':')
    plt.subplots_adjust(left=0.17, right=0.9, top=0.97, bottom=0.07, hspace=0.3)
    
    cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.2, aspect=40)
    cb1.set_label('Angular momentum gradient')
    
    plt.savefig(plot_dir+'ang_mom_grad_' + run + '.pdf', format='pdf')
    plt.close()
    

if __name__ == "__main__":
    
    ang_mom_grad_plot('sn_1_sst_zs')
    ang_mom_grad_plot('sn_1.000_zs_sst')


    
    
    
    