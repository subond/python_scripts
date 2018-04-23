# Locate maximum in zonal mean stretching term

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

def vort_eq_hm(run, lev=150, rotfac=1.0, period_fac=1.):
    
    rcParams['figure.figsize'] = 15, 6.25
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/crit_lat_test/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Load in vorticity budget term means
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_'+run+'.nc')
    
    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    omega = 7.2921150e-5 * rotfac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    
    div = gr.ddx(data.ucomp) + gr.ddy(data.vcomp)
    stretching_mean = -86400.**2. * vor * div
    stretching_hm = stretching_mean.mean('lon')
    
    mn_dic = month_dic(1)
    tickspace = np.arange(13,72,18) * period_fac
    labels = [mn_dic[(k+5)/6 ] for k in range(13, 72, 18)]
    levels = np.arange(-3.,3.1,0.25)
        
    stretching_max = np.unravel_index( stretching_hm.sel(pfull=lev).isel(lat=range(32,64)).argmin(), (data.pentad.size,32) )
    print stretching_max
    print data.lat[stretching_max[1] + 32]

    stretching_hm.sel(pfull=lev).plot.contourf(x='pentad', y='lat', levels=levels, extend = 'both', add_labels=False)
    #ax4.contour(data.xofyear, data.lat, abs_vort.T, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    plt.plot(data.pentad[stretching_max[0]], data.lat[stretching_max[1]+32], 'kx')
    plt.ylabel('Latitude')
    plt.xlabel('')
    plt.ylim(-60,60)
    plt.yticks(np.arange(-60,61,30))
    plt.xticks(tickspace,labels,rotation=25)
    plt.title('Vortex stretching', fontsize=17)
    plt.grid(True,linestyle=':')
    plt.tight_layout()  
    
    figname = 'vort_stretching_' + run + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
      
#vort_eq_hm('sn_1.000')
#vort_eq_hm('sn_2.000', period_fac=2.)
#vort_eq_hm('sn_0.500', period_fac=0.5)
vort_eq_hm('rt_2.000', rotfac=2.)
#vort_eq_hm('rt_0.500', rotfac=0.5)