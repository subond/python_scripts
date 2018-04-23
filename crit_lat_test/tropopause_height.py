# Estimate tropopause height using lapse rate method

"""
Load in the vorticity budget terms and produce a lat-time plot at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams

    
def lapse_rate(run, pentad, period_fac=1.0):
    
    #Load in data
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    g = 9.8
    Ra = 287.04
    
    mn_dic = month_dic(1)
    tickspace = np.arange(13,72,18) * period_fac
    labels = [mn_dic[(k+5)/6 ] for k in range(13, 72, 18)]
    levels = np.arange(-1.5,1.6,0.25)
    
    dlnTdp = gr.ddp(np.log(data.temp.mean('lon')))
    dTdz = data.pfull*100.*g/Ra * dlnTdp * 1000.

    #dTdz.mean('xofyear').plot.contourf(x='lat', y='pfull', extend = 'both', add_labels=False, yincrease=False, levels=np.arange(-10.,10.5,2.))
    #plt.xlabel('Latitude')
    #plt.ylabel('Pressure, hPa')
    #plt.grid(True,linestyle=':')
    #plt.tight_layout()  
    
    #figname = 'dTdz_' + run + '.pdf'
    #plt.savefig(plot_dir + figname, format='pdf')
    #plt.close()
    
    #ind = abs(dTdz.isel(pfull=range(10,20)) -2.).mean('xofyear').argmin('pfull')
    ind = abs(dTdz.isel(pfull=range(10,20)) -2.).isel(xofyear=pentad).argmin('pfull')
    
    trop_height = np.zeros([64,1])
    
    for i in range(0,64):
        #trop_height[i] = data.height.mean(('lon','xofyear'))[ind.values[i]+10,i]
        trop_height[i] = data.height.mean(('lon')).isel(xofyear=pentad)[ind.values[i]+10,i]
    
    print run, trop_height[32]
    
    return trop_height

plot_dir = '/scratch/rg419/plots/crit_lat_test/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

trop_height_sn10 = lapse_rate('sn_1.000', pentad=46)
trop_height_sn20 = lapse_rate('sn_2.000', pentad=85)
trop_height_sn05 = lapse_rate('sn_0.500', pentad=26)
trop_height_rt20 = lapse_rate('rt_2.000', pentad=59)
trop_height_rt05 = lapse_rate('rt_0.500', pentad=46)

plt.plot(trop_height_sn10)
plt.plot(trop_height_sn20)
plt.plot(trop_height_sn05)
plt.plot(trop_height_rt20)
plt.plot(trop_height_rt05)
plt.savefig(plot_dir + 'trop_heights.pdf', format='pdf')
plt.close()
