"""
Compare precip centroid and 0 mass streamfunction latitude
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#import sh
from data_handling import cell_area
#from pylab import rcParams
from physics import precip_centroid, model_constants as mc, mass_streamfunction, get_edge_psi
import scipy.interpolate as spint
from pylab import rcParams

month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 20
rcParams['text.usetex'] = True

def pcent_vs_0psi(run, period_fac=1., lonin=[-1.,361.]):
    
    #Load in data, add area to dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    psi_edge = get_edge_psi(run, thresh=0., period_fac=period_fac)
    
    p_cent = precip_centroid(run, lat_bound = 45., lonin=lonin)

    plt.plot(p_cent.sel(xofyear=psi_edge[0].xofyear), psi_edge[0], 'xk', alpha=0.7)
    plt.xlabel('Precipitation centroid')
    plt.ylabel('Latitude of 0 mass streamfunction')
    plt.grid(True,linestyle=':')
    plt.ylim([-40,40])
    plt.xlim([-40,40])
    plt.savefig('/scratch/rg419/plots/aht_work/pcent_0psi_'+run+'.pdf', format='pdf')
    plt.close()
    



pcent_vs_0psi('ap_2')
pcent_vs_0psi('ap_20')
pcent_vs_0psi('full_qflux', lonin=[60.,150.])
pcent_vs_0psi('sn_1.000')
pcent_vs_0psi('sn_2.000')
pcent_vs_0psi('sn_0.500')
pcent_vs_0psi('rt_0.500', period_fac=0.5)
pcent_vs_0psi('rt_2.000', period_fac=2.0)