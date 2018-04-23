"""
Function to plot meridional overturning strength vs northward extent

"""
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import numpy as np
import xarray as xr
from data_handling import time_means
from physics import mass_streamfunction
import sh
from finite_difference import cfd
from pylab import rcParams
import scipy.interpolate as spint

#import statsmodels.api as sm


    
rcParams['figure.figsize'] = 9, 6
rcParams['font.size'] = 20
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/rotation/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

lev = 500.
thresh = 60.

def get_edge_psi(run, lonin=[-1.,361.], sanity_check=True, lev=lev):
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180.)
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    #lats = data.lat[np.abs(data.lat) <= 35.]
    
    #Calculate mass streamfunction
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50., lons=lons)
    psi /= 1.e9
    
    f = spint.interp1d(psi.lat, psi, axis=0, fill_value='extrapolate')
    lats_new = np.arange(-90, 90.1, 0.5)
    psi_new = f(lats_new)
    psi_new = xr.DataArray(psi_new, coords=[lats_new, psi.xofyear, psi.pfull], dims=['lat', 'xofyear', 'pfull'])
    
    lats = psi_new.lat[np.abs(psi_new.lat) <= 35.]
    
    
    #Create array of zeros in lat-time shape
    edge_loc = np.zeros(psi_new.sel(pfull=lev).values.shape)
    # Mask where overturning magnitude is less than -50.
    edge_loc[psi_new.sel(pfull=lev).values <= -1.*thresh] = 1.
    edge_loc = xr.DataArray( edge_loc, psi_new.sel(pfull=lev).coords) 
    
    #Take diff along lat axis between 35 N and S and find location of first maximum in range
    eli = edge_loc.sel(lat=lats).diff('lat').argmin('lat')
    
    edge_loc_lat = psi_new.lat.sel(lat=lats)[eli]
    psi_min = -1.*psi_new.sel(pfull=lev).sel(lat=lats).min('lat')
    
    if sanity_check:
        #Plot as sanity check
        psi.sel(pfull=lev).plot.contourf(levels=np.arange(-300.,301.,50.))
        plt.plot(edge_loc_lat,'r')
        
        plt.figure(2)
        plt.plot(psi_min)
        
        plt.figure(3)
        plt.plot(edge_loc_lat[21:36],psi_min[21:36],'b')
        plt.plot(edge_loc_lat[50:70],psi_min[50:70],'g')
        plt.xlim([0,30])
        plt.show()
    
    return edge_loc_lat, psi_min


edge_ll_10, psi_min_10 = get_edge_psi('sn_1.000')
edge_ll_20, psi_min_20 = get_edge_psi('rt_2.000')
edge_ll_05, psi_min_05 = get_edge_psi('rt_0.500')



plt.plot(edge_ll_20, psi_min_20, 'bx', ms=10, mew=2)
plt.plot(edge_ll_05, psi_min_05, 'rx', ms=10, mew=2)
plt.plot(edge_ll_10, psi_min_10, 'kx', ms=10, mew=2)
ax = plt.gca()
#plt.xlabel('ITCZ latitude')
#plt.ylabel('Max 500 hPa $\Psi$, 10$^9$ kg/s')
plt.xlabel('$\phi_{0}$')
plt.ylabel('$\Psi_{max}$, 10$^9$ kg/s')
plt.xscale('log', subsx=[2,3,4,5])
plt.yscale('log', subsx=[2,3,4,5])
plt.xlim(1,40)
plt.ylim(30,600)
plt.xticks([1,5,10,20,40])
plt.yticks([30,75,150,300,600])
ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
ax.get_yaxis().set_major_formatter(tk.ScalarFormatter())
plt.tight_layout()  
figname = 'psi_latvsstrength_' + str(int(lev)) + '_' + str(int(thresh)) + '.pdf'
#figname = 'psi_latvsstrength.pdf'
plt.savefig(plot_dir + figname, format='pdf')
plt.close()
    
