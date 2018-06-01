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

    
def rossby_plot(run, ax, rot_fac=1., lev=150.):
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
    vor = (v_dx - u_dy).sel(pfull=lev)

    vor = make_sym(vor, asym=True)
    rossby = (-1.*vor/f).mean('lon')
    
    f1 = rossby.plot.contourf(ax=ax, x='xofyear', y='lat', levels=np.arange(-1.5,1.6,0.3), add_colorbar=False, add_labels=False)
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax, color='k', linewidth=2)
    ax.set_ylim([-60,60])
    ax.set_ylabel('Latitude')
    ax.set_xticks([12,24,36,48,60,72])
    ax.set_yticks([-60,-30,0,30,60])
    ax.set_xlabel('Pentad')
    ax.grid(True,linestyle=':')
    
    #plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.05)
    cb1=fig.colorbar(f1, ax=ax, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.2, aspect=40)
    cb1.set_label('Rossby number')


if __name__ == "__main__":
    
    # Set plotting directory
    plot_dir = '/scratch/rg419/plots/paper_2_figs/rossby/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)

    # Set figure parameters
    rcParams['figure.figsize'] = 6, 5
    rcParams['font.size'] = 14
    
    i=0
    rots=[0.5,0.75,1.,1.25,1.5,1.75,2.]
    for run in ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        rossby_plot(run, ax1, rot_fac = rots[i])
        plt.savefig(plot_dir+'rossby_no_' + run + '.pdf', format='pdf')
        plt.close()
        i=i+1
    
    
    