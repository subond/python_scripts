# Track the precipitation maximum over the year

from data_handling import month_dic
from physics import gradients as gr
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh


#rcParams['figure.figsize'] = 6, 10
rcParams['font.size'] = 20
rcParams['text.usetex'] = True
    

plot_dir = '/scratch/rg419/plots/seasons/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

def precip_max(run):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['totp'] = ((data.precipitation)*86400.).mean('lon')
    #lat_max = data.lat[data.totp.argmax('lat')].values
    
    lat_max = xr.DataArray(data.lat[data.totp.argmax('lat')].values, coords=[data.xofyear])
    
    ddt_lat = gr.ddt(lat_max)*86400.
        
    return lat_max, ddt_lat


lat_max_05, ddt_lat_05 = precip_max('sn_0.500')
lat_max_10, ddt_lat_10 = precip_max('sn_1.000')
lat_max_20, ddt_lat_20 = precip_max('sn_2.000')



mn_dic = month_dic(1)
tickspace = range(13,72,18)
labels = [mn_dic[(k+5)/6 ] for k in tickspace]    

plt.plot(lat_max_05.xofyear*2., lat_max_05, 'b')
plt.plot(lat_max_10.xofyear, lat_max_10, 'k')
plt.plot(lat_max_20.xofyear/2., lat_max_20, 'r')
plt.plot([40,40], [-30,30], '--k')
plt.xlim([0,73])
plt.xlabel('')
plt.ylabel('Latitude')
plt.xticks(tickspace, labels, rotation=25)
plt.grid(True,linestyle=':')

a = plt.axes([.18, .57, .3, .3])
plt.plot(lat_max_05.xofyear, lat_max_05, 'b')
plt.plot(lat_max_05.xofyear, lat_max_10[20:56], 'k')
plt.plot(lat_max_05.xofyear, lat_max_20[60:96], 'r')
plt.plot([20,20], [-30,30], '--k')
plt.xlim([0,37])
plt.xticks([])
plt.yticks([])
plt.ylabel('Latitude', fontsize=12)
plt.xlabel('Time', fontsize=12)

plt.savefig(plot_dir + 'precip_max_lat.pdf', format = 'pdf')
plt.close()