''' 
13/11/2017
Find precip centroid and overplot for different runs
'''

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from climatology import precip_centroid



rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 18
rcParams['text.usetex'] = True

def load_pcent_find_zero(data, lat_bound=45, lonin=[-1.,361.]):
    precip_centroid(data)
    
    pcent_sign = np.ones(data.p_cent.values.shape)
    
    pcent_sign[data.p_cent.values <= 0.] = 0.
    pcent_diff = np.diff(pcent_sign)
    loc = np.argmax(pcent_diff)
    len1 = np.max(pcent_sign.shape) - loc
    
    pcent_repos = np.zeros(data.p_cent.values.shape)
    pcent_repos[:len1] = data.p_cent[loc:]
    pcent_repos[len1:] = data.p_cent[:loc]
    
    pcent_repos = xr.DataArray(pcent_repos, coords=[data.xofyear.values], dims=['xofyear'])

    return pcent_repos



def plot_groups(data, plotname, key_labels=None, xlimits=None):
    
    n = len(data)
    for i in range(n):
        pcent = load_pcent_find_zero(data[i])
        plt.plot(pcent)
    
    if not xlimits==None:
        plt.xlim(xlimits)
    
    if not key_labels==None:
        plt.legend(key_labels)
    
    plt.xlabel('Time, pentads')
    plt.ylabel('Precipitation centroid')
    
    plt.savefig('/scratch/rg419/plots/itcz_plots/' + plotname +'.pdf', format='pdf')
    plt.close()
    
    

data_sn10 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')
data_sn05 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_0.500.nc')
data_sn20 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_2.000.nc')
data = [data_sn10, data_sn05, data_sn20]
plot_groups(data, 'pcent_seasons', ['360 day year', '180 day year', '720 day year'], xlimits=[0.,36.])


data_rt05 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_0.500.nc')
data_rt20 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_2.000.nc')
data = [data_sn10, data_rt05, data_rt20]
plot_groups(data, 'pcent_rotation', ['1*$\Omega$', '0.5*$\Omega$', '2*$\Omega$'], xlimits=[0.,36.])


data_ap2= xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ap_2.nc')
data_ap20 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ap_20.nc')
data = [data_ap2, data_sn10, data_ap20]
plot_groups(data, 'pcent_mlds', ['2m', '10m', '20m'], xlimits=[0.,36.])