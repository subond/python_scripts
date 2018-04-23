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

    
rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 25
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/clean_diags/ap_2/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

def get_edge_psi(run, lonin=[-1.,361.]):
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    lats = data.lat[np.abs(data.lat) <= 30.]

    psi = mass_streamfunction(data, a=6376.0e3, lons=lons, dp_in=50.)
    psi /= 1.e9
    
    cp = 287.04/2*7
    L = 2.500e6
    #data['mse'] = (cp*data.temp + L*data.sphum).sel(lon=lons).mean('lon')
    data['mse'] = (data.convection_rain + data.condensation_rain).sel(lon=lons).mean('lon')
    
    #mse_max_loc = data.lat[data.mse.sel(pfull=850.).argmax('lat').values].values
    mse_max_loc = data.lat[data.mse.argmax('lat').values].values
    
    psi_min = -1.*psi.sel(pfull=500.).sel(lat=lats).min('lat')
    
    #psi.sel(pfull=500.).plot.contourf(x='xofyear', y='lat')
    #plt.figure(2)
    #plt.plot( psi_min)
    #plt.figure(3)
    #plt.plot(mse_max_loc)
    #plt.show()
    
    return mse_max_loc, psi_min


mse_max_ap, psi_min_ap = get_edge_psi('ap_2')
mse_max_20, psi_min_20 = get_edge_psi('ap_20')
mse_max_fq, psi_min_fq = get_edge_psi('full_qflux', lonin=[60.,150.])

plt.plot(mse_max_ap)
plt.figure(2)
plt.plot(psi_min_ap)
#plt.figure(3)
#plt.plot(mse_max_fq)
#plt.figure(4)
#plt.plot(psi_min_fq)

plt.figure(5)
plt.plot(mse_max_ap, psi_min_ap,'kx')
#plt.plot(mse_max_20, psi_min_20,'bx')
#plt.plot(mse_max_fq, psi_min_fq,'rx')
#plt.xlim(0,40)
plt.xscale('log', subsx=[2,3,4,5])
plt.show()

low_lat = (edge_ll_ap <= 9.) & (edge_ll_ap >= 0.)
high_lat = (edge_ll_ap > 9.)
        
edge_low_lat = edge_ll_ap[low_lat]
psi_min_low_lat = psi_min_ap[low_lat]
    
edge_high_lat = edge_ll_ap[high_lat]
psi_min_high_lat = psi_min_ap[high_lat]
    
A = np.array([ edge_low_lat**(1./5.), np.ones(edge_low_lat.shape) ])
consts_low_lat = np.linalg.lstsq(A.T,psi_min_low_lat)[0] # obtaining the parameters
print consts_low_lat
    
A = np.array([ edge_high_lat**(3./4.), np.ones(edge_high_lat.shape) ])
consts_high_lat = np.linalg.lstsq(A.T,psi_min_high_lat)[0] # obtaining the parameters
print consts_high_lat
    
line_low = consts_low_lat[0]*np.arange(1.,9.)**(1./5.) + consts_low_lat[1] # regression line
line_high = consts_high_lat[0]*np.arange(10.,31.)**(3./4.) + consts_high_lat[1] # regression line

plt.plot(np.arange(1.,9.), line_low, 'k-', linewidth=2)
plt.plot(np.arange(10.,31.), line_high, 'k-', linewidth=2)
    
plt.plot(edge_ll_ap, psi_min_ap, 'kx', markersize=10)
plt.plot(edge_ll_20, psi_min_20, 'bx', markersize=10)
#plt.plot(edge_ll_fq, psi_min_fq, 'rx', markersize=10)
ax = plt.gca()
plt.xlabel('Northern edge of SH Hadley cell')
plt.ylabel('Peak strength of SH Hadley circulation')
plt.xscale('log', subsx=[2,3,4,5])
plt.xlim(1,30)
plt.ylim(50,500)
plt.xticks([1,5,10,20,30])
ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
plt.tight_layout()  
#figname = 'psi_latvsstrength_both.pdf'
figname = 'psi_latvsstrength_aps.pdf'
#plt.savefig(plot_dir + figname, format='pdf')
plt.close()
    
