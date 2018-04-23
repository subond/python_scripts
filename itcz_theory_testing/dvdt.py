"""
Evaluate and plot dvdt

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams

def dvdt_hm(run, lev=150, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 12, 5.5
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/itcz_theory_testing/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    dvdt = gr.ddt( data.vcomp.sel(lon=lons).mean('lon').sel(pfull=lev) ) * 86400.

    levels = np.arange(-20.,20.1,2.)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    # Six subplots
    plt.set_cmap('RdBu_r')
    #First plot
    f1=dvdt.plot.contourf(x='xofyear', y='lat', extend='both', add_colorbar=False, add_labels=False)
    plt.ylabel('Latitude')
    plt.title('dvdt', fontsize=17)
    plt.ylim(-60,60)
    plt.grid(True,linestyle=':')
    plt.yticks(np.arange(-60.,61.,30.))
    
    #Colorbar
    cb1=plt.colorbar(f1,use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
    
    if lonin == [-1.,361.]:
        figname = 'dvdt_' +run+ '.pdf'
    else:
        figname = 'dvdt_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()

dvdt_hm('ap_2')
