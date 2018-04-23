"""
Function to plot meridional overturning and absolute vorticity

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means
from physics import mass_streamfunction
import sh
from finite_difference import cfd
from pylab import rcParams


    
def upsi(run, months, before_after, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 6, 9
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    #data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
        
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    
    psi = mass_streamfunction(data, a=6376.0e3, lons=lons, dp_in=50.)
    psi /= 1.e9
    
    psi_before = psi[::-1,before_after[0]:before_after[0]+5,::-1].mean('xofyear').transpose()
    psi_after = psi[::-1,before_after[1]:before_after[1]+5,::-1].mean('xofyear').transpose()
    
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    abs_vort_before = (data.vor[before_after[0]:before_after[0]+5,:,:,:].mean('xofyear') + f)*86400.
    abs_vort_after = (data.vor[before_after[1]:before_after[1]+5,:,:,:].mean('xofyear') + f)*86400.
    
    # Two subplots
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    plt.set_cmap('RdBu_r')
    #First plot
    f1 = ax1.contourf(data.lat, data.pfull, abs_vort_before.sel(lon=lons).mean('lon'), levels=np.arange(-14,15,2), extend = 'both')
    ax1.contour(data.lat, data.pfull, psi_before, colors='k', levels=np.arange(-300.,301.,60.), linewidths=2)
    ax1.invert_yaxis()
    ax1.set_ylabel('Pressure, hPa')
    ax1.set_xlim(-60,60)
    ax1.grid(True,linestyle=':')
    
    #Second plot
    f2 = ax2.contourf(data.lat, data.pfull, abs_vort_after.sel(lon=lons).mean('lon'), levels=np.arange(-14,15,2), extend = 'both')
    ax2.contour(data.lat, data.pfull, psi_after, colors='k', levels=np.arange(-300.,301.,60.), linewidths=2)
    ax2.grid(True,linestyle=':')
    ax2.set_xlim(-60,60)
    ax2.invert_yaxis()
    ax2.set_ylabel('Pressure, hPa')
    ax2.set_xlabel('Latitude')
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0., hspace=0.2)
    #Colorbar
    cb1=f.colorbar(f2, ax=[ax1,ax2], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30)
    cb1.set_label('Absolute vorticity, day$^{-1}$')
    
    if lonin == [-1.,361.]:
        figname = 'upsi_vor.pdf'
    else:
        figname = 'upsi_vor_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir+figname, format='pdf')
    plt.close()
    
    


upsi('ap_2', [121,361], [30,39])
upsi('full_qflux', [121,481], [18,39])
upsi('full_qflux', [121,481], [18,39], lonin = [60.,150.])
#upsi('flat_qflux', [121,481], [18,44], lonin = [60.,150.])


