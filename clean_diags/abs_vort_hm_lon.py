"""
Evaluate and plot momentum budget at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, month_dic
import sh
from physics import gradients as gr
from pylab import rcParams
from windspharm.xarray import VectorWind

def abs_vort_hm(run, lev=150, filename='plev_daily', timeav='pentad', period_fac=1.,latin=25.):
    
    rcParams['figure.figsize'] = 9, 9
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
    names = [name_temp % m for m in range( 397, 409)  ]
        
    #read data into xarray 
    data = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    
    uwnd = data.ucomp.sel(pfull=lev)
    vwnd = data.vcomp.sel(pfull=lev)

    # Create a VectorWind instance to handle the computation of streamfunction and
    # velocity potential.
    w = VectorWind(uwnd, vwnd)

    # Compute the streamfunction and velocity potential.
    data['vor'], data['div'] = w.vrtdiv()
    
    #data = xr.open_dataset('/scratch/rg419/Data_moist/'+run+'climatologies/'+run+'.nc')

    #Coriolis
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    
    lat_hm = data.lat[np.argmin(np.abs(data.lat - latin))]
    
    abs_vort = (f + data.vor).sel(lat=lat_hm)*86400.
    
    levels = np.arange(0.,14.1,2.)
    
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    # Plot
    f1=abs_vort.plot.contourf(x='lon', y='time', extend='both', levels=levels, add_colorbar=False, add_labels=False)
    plt.set_cmap('inferno_r')
    plt.ylabel('Time')
    #plt.yticks(tickspace, labels, rotation=25)
    plt.xlabel('Longitude')
    plt.xlim(60,150)
    plt.ylim(240+33*360,90+33*360)
    plt.grid(True,linestyle=':')
    plt.tight_layout()  
    #Colorbar
    cb1=plt.colorbar(f1, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
    cb1.set_label('$day^{-1}$')
    
    figname = 'abs_vort_lon_hm.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
        
        
abs_vort_hm('ap_2')
abs_vort_hm('full_qflux')



