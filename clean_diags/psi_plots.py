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


    
def upsi(run, months, before_after, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 8, 7
    rcParams['font.size'] = 25
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
    
    psi_before = psi[::-1,before_after[0]:before_after[0]+4,::-1].mean('xofyear').transpose()
    psi_after = psi[::-1,before_after[1]:before_after[1]+4,::-1].mean('xofyear').transpose()
    
    plt.contour(data.lat, data.pfull, psi_before, colors='k', levels=np.arange(-300.,301.,50.), linewidths=2)
    plt.gca().invert_yaxis()
    plt.ylabel('Pressure, hPa')
    plt.xlabel('Latitude')
    plt.xlim(-60,60)
    plt.grid(True,linestyle=':')
    plt.tight_layout()
    plt.savefig(plot_dir+'psi_below.pdf', format='pdf')    
    plt.close()
    


    
    
#upsi('ap_2', [121,361], [30,39])
upsi('full_qflux', [121,481], [18,39])
#upsi('flat_qflux', [121,481], [18,44], lonin = [60.,150.])


