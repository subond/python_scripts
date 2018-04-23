"""
Plot overturning at 500 hPa at the Equator vs the interpolated precipitation centroid
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#import sh
from data_handling import cell_area
#from pylab import rcParams
from physics import precip_centroid, model_constants as mc, mass_streamfunction
import scipy.interpolate as spint
from pylab import rcParams

month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 20
rcParams['text.usetex'] = True

def pcent_vs_psi_eq(run, period_fac=1., lonin=[-1.,361.], do_plot=False):
    
    #Load in data, add area to dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        
    psi = mass_streamfunction(data, mc.a, dp_in=50., lons=lons)
    psi /= 1.e9
    psi_eq = psi.sel(pfull=500.).isel(lat=[31,32]).mean('lat')
        
    psi_eq.coords['month'] = np.mod( psi_eq.xofyear -1, 72.*period_fac) //(6*period_fac) +1    
    psi_eq_month = psi_eq.groupby('month').mean(('xofyear'))
    
    p_cent = precip_centroid(run, lat_bound = 45., lonin=lonin)
    p_cent.coords['month'] = np.mod( p_cent.xofyear -1, 72.*period_fac) //(6*period_fac) +1    
    p_cent_month = p_cent.groupby('month').mean(('xofyear'))
    
    if do_plot:
        plt.plot(p_cent, psi_eq, 'xk', alpha=0.7)
        plt.plot(p_cent_month, psi_eq_month, 'xk', ms=7, mew=2)
        for i in range(0, len(month_list)):
            plt.text(p_cent_month[i]+0.1, psi_eq_month[i]+0.1, month_list[i], fontsize=14)
        plt.xlabel('Precipitation centroid')
        plt.ylabel('500 hPa Mass Streamfunction at the Equator')
        plt.grid(True,linestyle=':')
        plt.ylim([-600,600])
        plt.xlim([-30,30])
        plt.savefig('/scratch/rg419/plots/aht_work/pcent_psi_'+run+'.pdf', format='pdf')
        plt.close()
    
    return p_cent, p_cent_month, psi_eq, psi_eq_month


def plot_pcent_psi(vars, varcolor='k', printmonths=False):
    #plt.plot(vars[0], vars[2], 'x', color=varcolor, alpha=0.7)
    plt.plot(vars[1], vars[3], 'x', color=varcolor, ms=7, mew=2)
    if printmonths:
        for i in range(0, len(month_list)):
            plt.text(vars[1][i]+0.1, vars[3][i]+0.1, month_list[i], fontsize=14, color = varcolor)


def plot3_pcent_psi(vars1, vars2, vars3, plotname):
    """Plot 3 sets of pcent vs psi on the same axis for comparison"""
    
    plot_pcent_psi(vars1, printmonths=True)
    plot_pcent_psi(vars2, varcolor='b')
    plot_pcent_psi(vars3, varcolor='r')
    plt.xlabel('Precipitation centroid')
    plt.ylabel('500 hPa Mass Streamfunction at the Equator')
    plt.grid(True,linestyle=':')
    plt.ylim([-600,600])
    plt.xlim([-30,30])
    plt.savefig('/scratch/rg419/plots/aht_work/pcentpsi_' + plotname +'.pdf', format='pdf')
    plt.close()
    

   
vars_ap2 = pcent_vs_psi_eq('ap_2')
vars_ap20 = pcent_vs_psi_eq('ap_20')
vars_full = pcent_vs_psi_eq('full_qflux', lonin=[60.,150.])
plot3_pcent_psi(vars_ap2, vars_ap20, vars_full, 'ap_to_full')

vars_sn10 = pcent_vs_psi_eq('sn_1.000')
vars_rt05 = pcent_vs_psi_eq('rt_0.500')
vars_rt20 = pcent_vs_psi_eq('rt_2.000')
plot3_pcent_psi(vars_sn10, vars_rt05, vars_rt20, 'rotation')

vars_sn05 = pcent_vs_psi_eq('sn_0.500', period_fac=0.5)
vars_sn20 = pcent_vs_psi_eq('sn_2.000', period_fac=2.)
plot3_pcent_psi(vars_sn10, vars_sn05, vars_sn20, 'seasons')

plot3_pcent_psi(vars_ap2, vars_sn10, vars_ap20, 'mld')


