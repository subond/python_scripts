"""
Reproduce moist versions of figures 8 and 9 from SB08 using data from the updated ap_2 run. 
Use Emanuel 1995 (On Thermally Direct Circulations...) angular momentum conserving equivalent 
potential temperature (12/11/2018)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, model_constants as mc, make_sym
from climatology import precip_centroid
from pylab import rcParams
from hadley_cell import mass_streamfunction, get_edge_psi


def find_change_lat(div_vt, sanity_check=False):
    # Find latitudes at which div(vt) changes from positive to negative, and set other values to NaN
    sign_changes = ((div_vt >=0.).astype(float).diff('lat') <0.).astype('float')*div_vt.lat
    sign_changes.values[sign_changes.values==0.] = np.nan
    
    # Find the location of the minimum div(vt), the sign change we want will be just south of this
    div_vt_min_loc = div_vt.lat.values[div_vt.argmin('lat').values]
    div_vt_min_loc = xr.DataArray(div_vt_min_loc, coords=[div_vt.xofyear], dims=['xofyear'])
    
    # Subtract the location of the sign changes from that of the minimum and exclude values north of the min
    distance_to_min = (div_vt_min_loc - sign_changes)
    distance_to_min.values[distance_to_min.values<0.] = np.nan
    
    # The closest value to the minimum is the one we want
    change_index = distance_to_min.argmin('lat')
    change_lat = xr.DataArray(div_vt.lat.values[change_index.values], coords=[div_vt.xofyear], dims=['xofyear'])
    
    if sanity_check==True:
        div_vt.plot.contourf(x='xofyear', y='lat')
        change_lat.plot()
        plt.show()
    return change_index, change_lat


def gross_stability(run, moist=False, i=1):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    convTtotheta=(1000./data.pfull)**mc.kappa
    
    theta = data.temp * convTtotheta
    theta_equiv = (data.temp + mc.L/mc.cp_air * data.sphum/(1-data.sphum)) * convTtotheta
    
    dthetady = gr.ddy(theta.mean('lon'), vector=False)
    dthetadp = gr.ddp(theta.mean('lon'))
    dthetady_equiv = gr.ddy(theta_equiv.mean('lon'), vector=False)
    dthetadp_equiv = gr.ddp(theta_equiv.mean('lon'))
    
    vdthetady_mean = data.vcomp.mean('lon') * dthetady
    wdthetadp_mean = data.omega.mean('lon') * dthetadp
    vdthetady_mean_equiv = data.vcomp.mean('lon') * dthetady_equiv
    wdthetadp_mean_equiv = data.omega.mean('lon') * dthetadp_equiv
    
    def column_int(var_in):
        var_int = mc.cp_air * var_in.sum('pfull')*5000./mc.grav
        return var_int
    
    div_vt_mean_int = -1. * column_int(wdthetadp_mean + vdthetady_mean)
    div_vt_mean_int_equiv = -1. * column_int(wdthetadp_mean_equiv + vdthetady_mean_equiv)
    
    vt_mean_int = column_int(data.vcomp.mean('lon') * theta.mean('lon'))
    vt_mean_int_equiv = column_int(data.vcomp.mean('lon') * theta_equiv.mean('lon'))
    #vt_mean_int.plot.contourf(x='xofyear', y='lat')
    
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    #psi /= 1.e9
    psi = np.abs(psi).max('pfull')
    #plt.figure(2)
    #psi.plot.contourf(x='xofyear', y='lat')
    
    gross_stab = (2.*np.pi * mc.a * np.cos(data.lat*np.pi/180.) * np.abs(vt_mean_int))/psi
    gross_moist_stab = (2.*np.pi * mc.a * np.cos(data.lat*np.pi/180.) * np.abs(vt_mean_int_equiv))/psi
    plt.figure(i)
    gross_moist_stab.plot.contourf(x='xofyear', y='lat', levels=np.arange(0., 2.e5, 2.e4))
    #plt.show()
    
    
gross_stability('sn_1.000_evap_fluxes_heattrans')
gross_stability('rt_0.500_heatbudg', i=2)
gross_stability('rt_0.750_heatbudg', i=3)
gross_stability('rt_1.250_heatbudg', i=4)
gross_stability('rt_1.500_heatbudg', i=5)
plt.show()