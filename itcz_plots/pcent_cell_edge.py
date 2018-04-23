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
    edge_loc, psi_max, psi_max_loc = get_edge_psi(data, thresh=0.)
    
    plt.plot(data.p_cent, edge_loc, symbol)
    
    
    

data_sn10 = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_0.500.nc')


data_zs = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/zs_sst.nc')


data_list = [xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_90.000.nc'),
             xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_91.000.nc'),
             xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_92.000.nc'),
             xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_93.000.nc'),
             xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_94.000.nc'),
             xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_95.000.nc'),
             xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_100.000.nc'),
             xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_105.000.nc'),
             xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_110.000.nc'),
             xr.open_dataset('/scratch/rg419/Data_moist/climatologies/ss_115.000.nc')
         ]


for data in data_list:
    pcent_vs_cell_edge(data)

pcent_vs_cell_edge(data_sn10, symbol='+r')
pcent_vs_cell_edge(data_zs, symbol='+g')

plt.show()