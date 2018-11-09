"""
Functions to look at absolute vorticity evolution in relation to precip centroid latitude (29/05/2018)
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

Rd = 287.04
Rv = 461.50
cp = Rd/2*7

def vdtdy_plot(run, ax, lev=850.):

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
    
    rho = lev*100./Rd/data.temp.sel(pfull=lev) / (1 + (Rv - Rd)/Rd*data.sphum.sel(pfull=lev))
    
    expansion_term = (1./rho/cp * data.omega.sel(pfull=lev)).mean('lon') * 86400.
    #expansion_term = (1./rho/cp ).mean('lon') * 86400.
    
    v = data.vcomp.sel(pfull=lev) # v
    T_dy = -86400. * gr.ddy( data.temp.sel(pfull=lev), vector=False)  # dTdy
    vdTdy = v.mean('lon') * T_dy.mean('lon') # [v][dTdy]
    
    w = data.omega.sel(pfull=lev)# w
    T_dp = -86400. * (gr.ddp(data.temp)).sel(pfull=lev)  # dTdp
    wdTdp = w.mean('lon') * T_dp.mean('lon') # [w][dTdp]
    
    dTdt = gr.ddt(data.temp).sel(pfull=lev).mean('lon')*86400.
    #dvtdy = -1.*(gr.ddy(data.vcomp * data.temp)).mean('lon').sel(pfull=lev)*86400.
    #dwtdp = -1.*(gr.ddp(data.omega * data.temp)).mean('lon').sel(pfull=lev)*86400.
    
    #dvtdy_colav = -1.*gr.ddy((data.vcomp * theta).mean('lon').sel(pfull=np.arange(50.,501.,50.)).mean('pfull')*5000./9.8)*1004.64
    
    #dvHdy_colav = -1.*gr.ddy((data.vcomp * (data.temp*1004.64 + data.sphum*2.500e6)).mean('lon').mean('pfull')*5000./9.8)
    
    
    #f1 = (dvtdy_colav).plot.contourf(ax=ax, x='xofyear', y='lat', add_labels=False, levels=np.arange(-30.,33.,3.), extend='both', add_colorbar=False)
    #f1 = (vdTdy).plot.contourf(ax=ax, x='xofyear', y='lat', add_labels=False, levels=np.arange(-5.,5.5,0.5), extend='both', add_colorbar=False)
    f1 = (vdTdy+wdTdp).plot.contourf(ax=ax, x='xofyear', y='lat', add_labels=False, levels=np.arange(-3.,3.1,0.5), extend='both', add_colorbar=False)
    #f1 = (dTdt).plot.contourf(ax=ax, x='xofyear', y='lat', add_labels=False, levels=np.arange(-0.2,0.21,0.02), extend='both', add_colorbar=False)
    #f1 = (vdTdy + wdTdp + expansion_term).plot.contourf(ax=ax, x='xofyear', y='lat', add_labels=False,  extend='both', add_colorbar=False)
    #f1 = (-1.*expansion_term-T_dp.mean('lon')).plot.contourf(ax=ax, x='xofyear', y='lat', add_labels=False,  extend='both', add_colorbar=False)
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax,color='k', linewidth=2)
    ax.set_xlabel('')
    ax.set_ylim([-60,60])
    ax.set_ylabel('Latitude')
    ax.set_xticks([12,24,36,48,60,72])
    ax.set_yticks([-60,-30,0,30,60])
    ax.grid(True,linestyle=':')
    
    return f1

#for run in ['rt_0.500', 'rt_0.750', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']:
for run in ['sine_sst_10m']:
#for run in ['ap_2']:
    plot_dir = '/scratch/rg419/plots/paper_2_figs/adv_tend_rot/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)

    rcParams['figure.figsize'] = 5.5, 7.
    rcParams['font.size'] = 14
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)

    f1 = vdtdy_plot('sn_1.000', ax1)
    vdtdy_plot(run, ax2)
    ax2.set_xlabel('Pentad')
    
    ax1.text(-10, 60., 'a)')
    ax2.text(-10, 60., 'b)')
    
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.0)

    cb1=fig.colorbar(f1, ax=(ax1, ax2), use_gridspec=True, orientation = 'horizontal',fraction=0.1, pad=0.1, aspect=40)
    cb1.set_label('Advective temperature tendency, Kday$^{-1}$')

    plt.savefig(plot_dir+'adv_tend_fig_850_nonflux_' + run + '.pdf', format='pdf')
    plt.close()
        

