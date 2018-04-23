"""
Function to plot meridional overturning, zonal wind speed, and eddy flux convergences

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means
from physics import mass_streamfunction
import sh
from finite_difference import cfd
from pylab import rcParams


    
def ampsi(run, before_after, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.], pres_mask=False):
    
    rcParams['figure.figsize'] = 10, 15
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
        
    psi = mass_streamfunction(data, a=6376.0e3, lons=lons, dp_in=50.)
    psi /= 1.e9
    
    ang_mom = (data.ucomp.sel(lon=lons).mean('lon') + omega*a*coslat)*a*coslat
    ang_mom = ang_mom/a**2/omega*15. #scale following SB08 to give neater units for plotting

    am_before = ang_mom[before_after[0]:before_after[0]+4,:,:].mean('xofyear')
    am_after = ang_mom[before_after[1]:before_after[1]+4,:,:].mean('xofyear')
        
    
    psi_before = psi[::-1,before_after[0]:before_after[0]+4,::-1].mean('xofyear').transpose()
    psi_after = psi[::-1,before_after[1]:before_after[1]+4,::-1].mean('xofyear').transpose()
    
    if pres_mask:
        pmask = xr.DataArray(np.ones(am_before.shape), am_before.coords)    
        p_max = (data.ps.sel(lon=lons).min('lon')).mean('xofyear')
        for i in range(0,64):
            a=[j for j in range(len(data.pfull)) if data.pfull[j] >= p_max[i]/100.]
            pmask[a,i] = float('NaN')
        am_before = am_before*pmask
        am_after = am_after*pmask
        psi_before = psi_before*pmask
        psi_after = psi_after*pmask
    
    # Two subplots
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    plt.set_cmap('RdBu_r')
    #First plot
    f1 = ax1.contour(data.lat, data.pfull, am_before, colors='0.75', extend = 'both', levels=np.arange(0.,17.))
    ax1.contour(data.lat, data.pfull, psi_before, colors='k', levels=np.arange(-300.,301.,50.))
    ax1.invert_yaxis()
    ax1.set_ylabel('Pressure, hPa')
    ax1.set_xlim(-60,60)
    ax1.grid(True,linestyle=':')
    
    #Second plot
    f2 = ax2.contour(data.lat, data.pfull, am_after, colors='0.75', extend = 'both', levels=np.arange(0.,17.))
    ax2.contour(data.lat, data.pfull, psi_after, colors='k', levels=np.arange(-300.,301.,50.))
    ax2.invert_yaxis()
    ax2.grid(True,linestyle=':')
    ax2.set_xlim(-60,60)
    ax2.set_ylabel('Pressure, hPa')
    ax2.set_xlabel('Latitude')
    
    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0.05, hspace=0.1)
    
    if lonin == [-1.,361.]:
        figname = 'ampsi.pdf'
    else:
        figname = 'ampsi_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
    
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    

#ampsi('ap_2', [30,39])
#ampsi('full_qflux', [18,39])
ampsi('full_qflux', [18,39], lonin = [60.,150.])
#ampsi('flat_qflux', [18,44])
#ampsi('flat_qflux', [18,44], lonin = [60.,150.])


