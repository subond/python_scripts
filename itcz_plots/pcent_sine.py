''' 
21/11/2017
Plot precip centroid vs overturning edge
'''

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from climatology import precip_centroid
from hadley_cell import get_edge_psi



rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 18
#rcParams['text.usetex'] = True

def pcent_vs_cell_edge(data, symbol='xk'):
    precip_centroid(data)
    x = np.linspace(0., 360., 72.)
    sinewave = -1.*np.sin(x*np.pi/180.)
    coswave = 3.*np.cos(3.*x*np.pi/180.) + 20.*np.cos(3. + x*np.pi/180.)
    
    
    plt.plot(sinewave, data.p_cent, symbol)
    plt.plot(sinewave, coswave)
    
    

data_sn10 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')

data_zs = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/zs_sst.nc')


pcent_vs_cell_edge(data_sn10, symbol='+r')
pcent_vs_cell_edge(data_zs, symbol='+g')

plt.show()