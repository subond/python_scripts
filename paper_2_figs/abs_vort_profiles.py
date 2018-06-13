"""
Plot absolute vorticity profiles for the control and zs runs

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_mse_plot
from pylab import rcParams
from hadley_cell import mass_streamfunction
from data_handling_updates import make_sym, gradients as gr, model_constants as mc

plot_dir = '/scratch/rg419/plots/paper_2_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
rcParams['figure.figsize'] = 5.5, 7
rcParams['font.size'] = 14

#fig = plt.figure()
#ax1 = fig.add_subplot(111)    
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)


def vor_plot(run, ax, tf, linestyle='-', rot_fac=1., lev=150.):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    # Calculate vorticity
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy + f).sel(pfull=lev)*86400.
    
    div = (gr.ddx(data.ucomp) + gr.ddy(data.vcomp)).sel(pfull=lev)
    
    vor = (make_sym(vor, asym=True)).mean('lon')
    div = (make_sym(div)).mean('lon')*8640.
    
    m = ((mc.omega * mc.a**2. * np.cos(data.lat*np.pi/180.)**2. + data.ucomp.mean('lon') * mc.a * np.cos(data.lat*np.pi/180.)).sel(pfull=lev))/(mc.omega * mc.a**2.)

    #dvordt = gr.ddt(vor.mean('lon'))*86400.
    #stretching = (-86400. * vor * div).mean('lon')
    #adv = -86400. * (data.vcomp.sel(pfull=lev) * gr.ddy(vor, vector=False)).mean('lon')
    #dvordy = gr.ddy(vor, vector=False).mean('lon')*8640.
    
    #dvordy[tf[0]:tf[1],:].mean('xofyear').plot(ax=ax, color='r', linestyle=linestyle)
    #m[:,tf[0]:tf[1]].mean('xofyear').plot(ax=ax, color='r', linestyle=linestyle)
    vor[tf[0]:tf[1],:].mean('xofyear').plot(ax=ax, color='k', linestyle=linestyle)
    #stretching[tf[0]:tf[1],:].mean('xofyear').plot(ax=ax, color='b', linestyle=linestyle)
    #adv[tf[0]:tf[1],:].mean('xofyear').plot(ax=ax, color='r', linestyle=linestyle)
    #(stretching + adv)[tf[0]:tf[1],:].mean('xofyear').plot(ax=ax, color='r', linestyle=linestyle)
    
    ax.set_title('')
    ax.set_xlabel('')
    #ax.set_ylim(0.5,1.)
    ax.set_ylim(-10,10)
    ax.grid(True,linestyle=':')
    ax.set_ylabel('Absolute vorticity, day$^{-1}$')

vor_plot('sn_1.000', ax=ax1, tf=[31,35])
vor_plot('sn_1_sst_zs', ax=ax1, tf=[31,35], linestyle='--')
vor_plot('sn_1.000', ax=ax2, tf=[38,42])
vor_plot('sn_1_sst_zs', ax=ax2, tf=[38,42], linestyle='--')
vor_plot('sn_1.000', ax=ax3, tf=[45,49])
vor_plot('sn_1_sst_zs', ax=ax3, tf=[45,49], linestyle='--')

ax1.legend(['control','control-zs'], loc='upper left', fontsize=10)
ax3.set_xlim(-35,35)
ax3.set_xticks(np.arange(-30,31,15))
ax3.set_xlabel('Latitude')

ax1.text(-52, 10., 'a)')
ax2.text(-52, 10., 'b)')
ax3.text(-52, 10., 'c)')

ax1.set_title('Pentad 32-35', fontsize=14)
ax2.set_title('Pentad 39-42', fontsize=14)
ax3.set_title('Pentad 46-49', fontsize=14)

plt.subplots_adjust(left=0.20, right=0.95, top=0.95, bottom=0.07)

plt.savefig(plot_dir+'abs_vort_profiles.pdf', format='pdf')
plt.close()        
