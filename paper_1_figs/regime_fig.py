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

import statsmodels.api as sm


    
rcParams['figure.figsize'] = 9, 6
rcParams['font.size'] = 20
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/paper_1_figs/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

lev = 500.
thresh = 120.

def get_edge_psi(run, lonin=[-1.,361.], sanity_check=True, lev=lev):
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180)
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    lats = data.lat[np.abs(data.lat) <= 35.]
    
    #Calculate mass streamfunction
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50., lons=lons)
    psi /= 1.e9
    
    #Create array of zeros in lat-time shape
    edge_loc = np.zeros(psi.sel(pfull=lev).values.shape)
    # Mask where overturning magnitude is less than -50.
    edge_loc[psi.sel(pfull=lev).values <= -1.*thresh] = 1.
    edge_loc = xr.DataArray( edge_loc, psi.sel(pfull=lev).coords) 
    
    #Take diff along lat axis between 35 N and S and find location of first maximum in range
    eli = edge_loc.sel(lat=lats).diff('lat').argmin('lat')
    
    edge_loc_lat = data.lat.sel(lat=lats)[eli]
    psi_min = -1.*psi.sel(pfull=lev).sel(lat=lats).min('lat')
    
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


edge_ll_ap, psi_min_ap = get_edge_psi('ap_2')
edge_ll_20, psi_min_20 = get_edge_psi('ap_20')
edge_ll_fq, psi_min_fq = get_edge_psi('full_qflux', lonin=[60.,150.])

def fit_power_law(edge_ll, psi_min, maxmin):
    
    lat = (edge_ll <= maxmin[1]) & (edge_ll > maxmin[0])
    edge = np.log( edge_ll[lat] )
    psi_min = np.log( psi_min[lat] )
    
    A = np.array([ edge, np.ones(edge.shape) ])
    #consts = np.linalg.lstsq( A.T, psi_min )[0] # obtaining the parameters
    #resids = np.linalg.lstsq( A.T, psi_min )[1] # obtaining the parameters


    model = sm.OLS(psi_min, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    
    #print result.summary()
    
    print '=== Coeffs ==='
    print np.exp(consts[1]), consts[0]
    print '=== Std Errs ==='
    print 2*std_err[1]*np.exp(consts[1]), 2*std_err[0]
    
    line = np.exp(consts[1]) * np.arange(maxmin[0],maxmin[1]+1)**(consts[0])
    
    return line, consts

edge_ll_aps = np.zeros((144,))
edge_ll_aps[0:72] = edge_ll_ap
edge_ll_aps[72:144] = edge_ll_20

psi_min_aps = np.zeros((144,))
psi_min_aps[0:72] = psi_min_ap
psi_min_aps[72:144] = psi_min_20

line_low_aps, consts_low = fit_power_law(edge_ll_aps, psi_min_aps, maxmin=[0.,9.])
line_high_aps, consts_high = fit_power_law(edge_ll_aps, psi_min_aps, maxmin=[9.,30.])

line_low_fq, consts_low = fit_power_law(edge_ll_fq.values, psi_min_fq.values, maxmin=[0.,9.])
line_high_fq, consts_high = fit_power_law(edge_ll_fq.values, psi_min_fq.values, maxmin=[9.,30.])

plt.plot(np.arange(0.,10.), line_low_aps, 'k:', linewidth=2)
plt.plot(np.arange(9.,31.), line_high_aps, 'k:', linewidth=2)

plt.plot(np.arange(0.,10.), line_low_fq, 'r:', linewidth=2)
plt.plot(np.arange(9.,31.), line_high_fq, 'r:', linewidth=2)

plt.plot(edge_ll_20, psi_min_20, 'bx', ms=10, mew=2)
plt.plot(edge_ll_fq, psi_min_fq, 'rx', ms=10, mew=2)
plt.plot(edge_ll_ap, psi_min_ap, 'kx', ms=10, mew=2)
ax = plt.gca()
#plt.xlabel('ITCZ latitude')
#plt.ylabel('Max 500 hPa $\Psi$, 10$^9$ kg/s')
plt.xlabel('$\phi_{0}$')
plt.ylabel('$\Psi_{max}$, 10$^9$ kg/s')
plt.xscale('log', subsx=[2,3,4,5])
plt.yscale('log', subsx=[2,3,4,5])
plt.xlim(1,30)
plt.ylim(175,500)
plt.xticks([1,5,10,20,30])
plt.yticks([200,300,400,500])
ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
ax.get_yaxis().set_major_formatter(tk.ScalarFormatter())
plt.tight_layout()  
figname = 'psi_latvsstrength_' + str(int(lev)) + '_' + str(int(thresh)) + '.pdf'
#figname = 'psi_latvsstrength.pdf'
plt.savefig(plot_dir + figname, format='pdf')
plt.close()
    
