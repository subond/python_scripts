"""
Load in the vorticity budget terms and produce a lat-time plot at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from physics import gradients as gr
from pylab import rcParams

def vort_eq_hm(run, bef_aft, lev=150, lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 10, 6.25
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/paper_1_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Also load climatological data so that transient eddies can be calculated 
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_uv'+run+'.nc')
    print 'climatology loaded'
        
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    
    abs_vort = vor.sel(lon=lons).mean('lon')*86400.
    
    print 'starting plotting'
    
    np.abs(vor.sel(lon=lons, xofyear = range(bef_aft[0],bef_aft[0]+4)).mean(('lon','xofyear'))*86400.).plot(color='k')
    np.abs(vor.sel(lon=lons, xofyear = range(bef_aft[1],bef_aft[1]+4)).mean(('lon','xofyear'))*86400.).plot(color='r')
    plt.xlim([-45,45])
    plt.ylim([-12,12])
    plt.grid(True,linestyle=':')
    

    if lonin == [-1.,361.]:
        figname = 'vort_line_' + run + '.pdf'
    else:
        figname = 'vort_line_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'

    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
      
vort_eq_hm('ap_2', [31,40])
vort_eq_hm('full_qflux', [19,40], lonin=[60.,150.])


