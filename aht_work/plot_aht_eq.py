"""
Calculate AHT_EQ, the atmospheric heat transport across the equator. Calculation based on Donohoe et al. 2012 (J. Clim.)
Plot over the year

"""

import matplotlib.pyplot as plt
#import xarray as xr
from data_handling import month_dic
from physics import aht_eq


def plot_aht_eq(run):
    
    ahteq, swabs, olr, shf, stor = aht_eq(run)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    ahteq.plot(color='m')
    (-1.*stor).plot(color='c')
    (-1.*olr).plot(color='g')
    shf.plot(color='b')
    swabs.plot(color='r')
    #plt.xlim((1,72))
    #plt.xlabel('')
    #plt.xticks(tickspace,labels,rotation=25)
    plt.ylabel('Atmospheric heat transport at the equator, PW')
    plt.savefig('/scratch/rg419/plots/aht_work/aht_eq_'+run+'.pdf', format='pdf')
    plt.close()
    
    

#ahteq = plot_aht_eq('rt_2.000')
#ahteq = plot_aht_eq('rt_0.500')
#ahteq = plot_aht_eq('sn_1.000')
ahteq = plot_aht_eq('sn_2.000')
ahteq = plot_aht_eq('sn_0.500')
#ahteq = plot_aht_eq('full_qflux_altevap2')