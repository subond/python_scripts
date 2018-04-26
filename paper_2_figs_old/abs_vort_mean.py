"""
Take mean of absolute vorticity magnitude over the Hadley cell

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from hadley_cell import get_streamfun_zeros
from data_handling_updates import gradients as gr
from climatology import precip_centroid
from pylab import rcParams


def abs_vort_mean(run, rot_fac=1., lev=150.):

    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    precip_centroid(data)
    
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (v_dx - u_dy + f).sel(pfull=lev)
    
    # Take gradients of vorticity
    dvordx = gr.ddx(vor)
    dvordy = gr.ddy(vor, vector=False)
    
    # Horizontal material derivative
    horiz_md_mean = (-86400.**2. * (data.ucomp.sel(pfull=lev) * dvordx + data.vcomp.sel(pfull=lev) * dvordy)).mean('lon')
    
    # Calculate divergence and stretching term
    div = gr.ddx(data.ucomp.sel(pfull=lev)) + gr.ddy(data.vcomp.sel(pfull=lev))
    stretching_mean = (-86400.**2. * vor * div).mean('lon')
    
    vor = vor.mean('lon')
    vor_mag = np.abs(vor)
    stretching_mag = stretching_mean
    stretching_mag[:,0:32] = stretching_mag[:,0:32]*-1.
    
    lat_itcz, lat_nh, lat_sh = get_streamfun_zeros(data, lev=lev)
    
    def take_cell_average(var, lats_min, lats_max):
        cell_av = np.zeros(len(data.xofyear.values))
        for i in range(len(data.xofyear.values)):
            lats = [data.lat[j] for j in range(len(data.lat)) if data.lat[j] > lats_min[i] and data.lat[j] < lats_max[i]]        
            cell_av[i] = var[i].sel(lat=lats).mean('lat')
        return cell_av
    
    vor_av_sh = take_cell_average(vor_mag, lat_sh, lat_itcz)
    stretching_av_sh = take_cell_average(stretching_mag, lat_sh, lat_itcz)
    
    plt.plot(stretching_av_sh, 'k')
    plt.figure(2)
    plt.plot(data.p_cent, stretching_av_sh, 'kx')
    plt.show()
    
abs_vort_mean('rt_0.500', rot_fac=0.5)
abs_vort_mean('sn_1.000')
abs_vort_mean('rt_2.000', rot_fac=2.)