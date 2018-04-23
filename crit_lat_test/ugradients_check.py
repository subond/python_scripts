# Compare dv/dy with vmax/half cell width

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

def grads_hm(run, lev=150, rotfac=1.0, period_fac=1.):
    
    rcParams['figure.figsize'] = 15, 6.25
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/crit_lat_test/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Load in data
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_'+run+'.nc')
    
    # Evaluate gradient
    dvdy = gr.ddy(data.vcomp).mean('lon')
    
    mn_dic = month_dic(1)
    tickspace = np.arange(13,72,18) * period_fac
    labels = [mn_dic[(k+5)/6 ] for k in range(13, 72, 18)]
    levels = np.arange(-1.5,1.6,0.25)

    dvdy.sel(pfull=lev).plot.contourf(x='pentad', y='lat', extend = 'both', add_labels=False, levels=np.arange(-0.000010,0.0000105,0.000001))
    plt.ylabel('Latitude')
    plt.xlabel('')
    plt.ylim(-60,60)
    plt.yticks(np.arange(-60,61,30))
    plt.xticks(tickspace,labels,rotation=25)
    plt.title('dv/dy, /s', fontsize=17)
    plt.grid(True,linestyle=':')
    plt.tight_layout()  
    
    figname = 'dvdy_' + run + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    data.vcomp.mean('lon').sel(pfull=lev).plot.contourf(x='pentad', y='lat', extend = 'both', levels=np.arange(-10.,10.5,1.), add_labels=False)
    plt.ylabel('Latitude')
    plt.xlabel('')
    plt.ylim(-60,60)
    plt.yticks(np.arange(-60,61,30))
    plt.xticks(tickspace,labels,rotation=25)
    plt.title('v, m/s', fontsize=17)
    plt.grid(True,linestyle=':')
    plt.tight_layout()  
    
    figname = 'v_' + run + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    print run, data.vcomp.mean('lon').sel(pfull=lev).min().values
      
grads_hm('sn_1.000')
grads_hm('sn_2.000', period_fac=2.)
grads_hm('sn_0.500', period_fac=0.5)
grads_hm('rt_2.000', rotfac=2.)
grads_hm('rt_0.500', rotfac=0.5)