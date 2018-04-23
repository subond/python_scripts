"""
Function to plot meridional overturning, zonal wind speed, and eddy flux convergences

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, rolling_mean
from physics import mass_streamfunction, gradients as gr, model_constants as mc
import sh
from pylab import rcParams
import scipy.interpolate as spint

    
rcParams['figure.figsize'] = 12, 7
rcParams['font.size'] = 18
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/time_series/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def interp_field(field, lats = np.arange(-90, 90.1, 0.1)):

    f = spint.interp1d(field.lat, field, axis=1, fill_value='extrapolate', kind='cubic')
    field_int = f(lats)
    field_int = xr.DataArray(field_int, coords=[field.xofyear, lats], dims=['xofyear', 'lat'])
    
    return field_int    


def max_div(run, lats = np.arange(-90, 90.1, 0.1), rotfac=1.0):
    
    #Open data
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_'+run+'.nc')    
    
    #Calculate vorticity, divergence, and pentad mean vortex stretching    
    
    f = 2 * mc.omega * rotfac * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    
    div = gr.ddx(data.ucomp) + gr.ddy(data.vcomp)
    stretching_mean = -86400.**2. * vor * div
    
    # Select upper level and take zonal mean
    # Take a 20 day window rolling mean, and interpolate onto a 0.1 deg lat grid.
    # Find max values of divergence and stretching, and locate lat of max divergence
    div_hm = div.mean('lon').sel(pfull=150.)
    div_rm = rolling_mean(div_hm, 4, tcoord='pentad')
    div_int = interp_field(div_rm)
    div_max = div_int.max('lat') * 86400.
    div_argmax = div_int.lat.values[div_int.argmax('lat').values]
    div_argmax = xr.DataArray(div_argmax, coords=[div_int.xofyear], dims=['xofyear'])
    

    stretching_hm = stretching_mean.mean('lon').sel(pfull=150.)
    stretching_hm[:,32:] = stretching_hm[:,32:]*-1.
    stretching_rm = rolling_mean(stretching_hm, 4, tcoord='pentad')
    stretching_int = interp_field(stretching_rm)
    stretching_max = stretching_int.max('lat')
    
    # Plot up max divergence and vortex stretching tendency in top panel, latitude of max divergence below
    #fig, ax1 = plt.subplots()
    fig, (ax1, ax3) = plt.subplots(2, sharex=True)
    div_max.plot(ax=ax1, color='k')
    ax1.set_xlabel('')
    ax1.set_title('')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('Divergence, day$^{-1}$')
    ax1.tick_params('y')
    ax1.grid(True,linestyle=':')
    ax1.set_xlim(0,max(div_max.xofyear)+1)

    ax2 = ax1.twinx()
    stretching_max.plot(color='b', ax=ax2)
    ax2.set_title('')
    ax2.set_ylabel('Vortex stretching, day$^{-2}$', color='b')
    ax2.tick_params('y', colors='b')
    ax2.set_xlim(0,max(div_max.xofyear)+1)
    
    div_argmax.plot(ax=ax3, color='k')
    ax3.set_title('')
    ax3.set_xlabel('Pentad')
    ax3.set_ylabel('Lat of max divergence')
    ax3.grid(True,linestyle=':')
    ax3.set_xlim(0,max(div_max.xofyear)+1)

    fig.tight_layout()
    plt.savefig(plot_dir + 'div_tseries_' + run + '.pdf', format='pdf')
    plt.close()    
    

#max_div('sn_0.500')
#max_div('sn_1.000')
#max_div('sn_2.000')
max_div('rt_0.500', rotfac=0.5)
max_div('rt_2.000', rotfac=2.)

