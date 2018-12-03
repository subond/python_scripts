"""
Functions to look at absolute vorticity evolution in relation to precip centroid latitude (19/04/2018)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, make_sym, cell_area, model_constants as mc
from climatology import precip_centroid
from pylab import rcParams
import scipy.interpolate as spint
from hadley_cell import mass_streamfunction

    
def rossby_plot(run, ax, rot_fac=1., lev=200., type='vor_only', plottype='rossby'):
    '''Plot dvordt or 1/vor * dvordt'''
    
    #Load data
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    #Calculate psi to overplot
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    psi = make_sym(psi, asym=True)
    
    # Make precip symmetric and find the precip centroid
    data['precipitation'] = make_sym(data.precipitation)
    data['ucomp'] = make_sym(data.ucomp)
    data['vcomp'] = make_sym(data.vcomp, asym=True)
    data['omega'] = make_sym(data.omega)
    
    precip_centroid(data)
    
    # Calculate vorticity
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    u_dp = gr.ddp(data.ucomp)  # dudp
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy).sel(pfull=lev)
    metric = (data.ucomp/mc.a * np.tan(data.lat *np.pi/180)).sel(pfull=lev)
    vertical_term = (-1.*data.omega/data.vcomp * u_dp).sel(pfull=lev)
    #vor = make_sym(vor, asym=True)
    #metric = make_sym(metric, asym=True)
    #vertical_term = make_sym(vertical_term, asym=True)
    if type=='vor_only':
        rossby = (-1.*vor/f).mean('lon')
    elif type=='metric':
        rossby = (-1.*(metric+vor)/f).mean('lon')
    elif type=='vertical':
        rossby = (-1.*(vor + vertical_term)/f).mean('lon')
    elif type=='full':
        rossby = (-1.*(metric + vor + vertical_term)/f).mean('lon')
    levels=np.arange(-0.9,1.0,0.3)
    if plottype=='drodt':
        rossby = gr.ddt(rossby) * 84600.
        levels=np.arange(-0.1,0.1,0.01)
    f1 = rossby.plot.contourf(ax=ax, x='xofyear', y='lat', levels=levels, add_colorbar=False, add_labels=False, extend='both')
    psi.sel(pfull=lev).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=lev).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=lev).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax, color='k', linewidth=2)
    ax.set_ylim([-30,30])
    ax.set_ylabel('Latitude')
    ax.set_xticks([12,24,36,48,60,72])
    ax.set_yticks([-30,-15,0,15,30])
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
    rots=[0.5,0.75,1.,1.25,1.5,1.75,2.,1.,1.,1.,1., 0.5,0.75]
    for run in ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000', 'mld_2.5', 'mld_5', 'mld_15', 'mld_20', 'rt_0.500_5', 'rt_0.750_5']:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        rossby_plot(run, ax1, rot_fac = rots[i])
        plt.savefig(plot_dir+'rossby_no_' + run + '_vor_only.pdf', format='pdf')
        plt.close()
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        rossby_plot(run, ax1, rot_fac = rots[i], type='metric')
        plt.savefig(plot_dir+'rossby_no_' + run + '_metric.pdf', format='pdf')
        plt.close()
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        rossby_plot(run, ax1, rot_fac = rots[i], type='full')
        plt.savefig(plot_dir+'rossby_no_' + run + '_full.pdf', format='pdf')
        plt.close()
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        rossby_plot(run, ax1, rot_fac = rots[i], type='vertical')
        plt.savefig(plot_dir+'rossby_no_' + run + '_vertical.pdf', format='pdf')
        plt.close()
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        rossby_plot(run, ax1, rot_fac = rots[i], type='full', plottype='drodt')
        plt.savefig(plot_dir+'rossby_no_' + run + '_drodt.pdf', format='pdf')
        plt.close()
        
        i=i+1
    
    
    