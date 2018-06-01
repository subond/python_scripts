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


def vdtdy_plot(run, lev=850.):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    
    try:
        data['precipitation'] = make_sym(data.precipitation)
    except:
        data['precipitation'] = data.convection_rain + data.condensation_rain
        data['precipitation'] = make_sym(data.precipitation)
    
    precip_centroid(data)
    
    convTtotheta=(1000./data.pfull)**(2./7.)
    
    theta = data.temp*convTtotheta
    #rcParams['figure.figsize'] = 5.5, 4.3
    #theta.mean('lon').sel(xofyear=40).plot.contourf(x='lat', y='pfull', yincrease=False)
    #plt.show()
    dvtdy = -1.*(gr.ddy(data.vcomp * theta)).mean('lon').sel(pfull=lev)*86400.
    dwtdp = -1.*(gr.ddp(data.omega * theta)).mean('lon').sel(pfull=lev)*86400.
    
    adv_tend_grad = gr.ddy(dvtdy+dwtdp, vector=False)*10.**5.
    
    dTdt = (gr.ddt(theta)).mean('lon').sel(pfull=lev)*86400./theta.mean('lon').sel(pfull=lev)
    
    dvtdy_colav = -1.*gr.ddy((data.vcomp * theta).mean('lon').mean('pfull')*5000./9.8)*1004.64
    #dwtdp = -1.*gr.ddp((data.omega * theta).mean('lon').sum('pfull')*50./9.8)
    dTdt_colav = gr.ddt(theta.mean('lon').sum('pfull')*5000./9.8)*1004.64
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/dvtdy/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #rcParams['figure.figsize'] = 5.5, 4.3
    rcParams['figure.figsize'] = 5, 10.5
    rcParams['font.size'] = 14
    #rcParams['font.size'] = 16

    #fig = plt.figure()
    #ax1 = fig.add_subplot(111)    
    fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1)
    
    f1 = (dvtdy+dwtdp).plot.contourf(ax=ax1, x='xofyear', y='lat', add_labels=False, levels=np.arange(-50.,55.,5.), extend='both')
    psi.sel(pfull=500).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax1, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax1,color='k', linewidth=2)
    ax1.set_ylim([-60,60])
    ax1.set_ylabel('Latitude')
    ax1.set_xticks([12,24,36,48,60,72])
    ax1.set_yticks([-60,-30,0,30,60])
    ax1.set_xlabel('Pentad')
    ax1.grid(True,linestyle=':')
    
    f1 = adv_tend_grad.plot.contourf(ax=ax2, x='xofyear', y='lat', add_labels=False, levels=np.arange(-10.,11.,1.),  extend='both')
    psi.sel(pfull=500).plot.contour(ax=ax2, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax2, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax2, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax2,color='k', linewidth=2)
    ax2.set_ylim([-60,60])
    ax2.set_ylabel('Latitude')
    ax2.set_xticks([12,24,36,48,60,72])
    ax2.set_yticks([-60,-30,0,30,60])
    ax2.set_xlabel('Pentad')
    ax2.grid(True,linestyle=':')
    
    #f1 = dvtdy_colav.plot.contourf(ax=ax3, x='xofyear', y='lat', add_labels=False, levels=np.arange(-50.,51.,5.), extend='both')
    #psi.sel(pfull=500).plot.contour(ax=ax3, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    #psi.sel(pfull=500).plot.contour(ax=ax3, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    #psi.sel(pfull=500).plot.contour(ax=ax3, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    #data.p_cent.plot.line(ax=ax3,color='k', linewidth=2)
    #ax3.set_ylim([-60,60])
    #ax3.set_ylabel('Latitude')
    #ax3.set_xticks([12,24,36,48,60,72])
    #ax3.set_yticks([-60,-30,0,30,60])
    #ax3.set_xlabel('Pentad')
    #ax3.grid(True,linestyle=':')
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.05)

   # cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.2, aspect=40)
   # cb1.set_label('Advective temperature tendency, Kday$^{-1}$')

    plt.savefig(plot_dir+'adv_tend_' + run + '.pdf', format='pdf')
    plt.close()
        


for run in ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']: 
    vdtdy_plot(run)

for run in ['mld_2.5','mld_5','mld_15','mld_20']: 
    print(run)
    vdtdy_plot(run)

for run in ['sn_0.500','sn_2.000','sn_3.000','sn_4.000']: 
    print(run)
    vdtdy_plot(run)
