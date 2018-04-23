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


    
def psi_plot(run, period_fac=1.,lonin=[-1.,361.], sanity_check=False):
    
    rcParams['figure.figsize'] = 10, 7
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
    
    edge_loc = np.zeros(psi.sel(pfull=500).values.shape)
    edge_loc[psi.sel(pfull=500).values <= -50.] = 1.
    edge_loc = xr.DataArray( edge_loc, psi.sel(pfull=500).coords) 
        
    eli = edge_loc.diff('lat')[::-1,:].argmax('lat')
    
    edge_loc_lat = data.lat[eli]
    edge_ll_summer = edge_loc_lat[abs(data.lat[eli]) < 30.]
    
    psi_min = psi.sel(pfull=500.).min('lat')
    psi_min_summer = -1.*psi_min[abs(data.lat[eli]) < 30.]
    
    if sanity_check:
        psi.sel(pfull=500).plot.contourf(levels=np.arange(-300.,301.,50.))
        plt.plot(edge_loc_lat)
        plt.figure(2)
        plt.plot(psi_min)
        plt.plot(psi_min_summer,'r')
        plt.show()
    

    low_lat = (edge_ll_summer <= 9.) & (edge_ll_summer >= 0.)
    high_lat = (edge_ll_summer > 9.)
        
    edge_low_lat = edge_ll_summer[low_lat]
    psi_min_low_lat = psi_min_summer[low_lat]
    
    edge_high_lat = edge_ll_summer[high_lat]
    psi_min_high_lat = psi_min_summer[high_lat]
    
    A = np.array([ edge_low_lat**(1./5.), np.ones(edge_low_lat.shape) ])
    consts_low_lat = np.linalg.lstsq(A.T,psi_min_low_lat)[0] # obtaining the parameters
    print consts_low_lat
    
    A = np.array([ edge_high_lat**(3./4.), np.ones(edge_high_lat.shape) ])
    consts_high_lat = np.linalg.lstsq(A.T,psi_min_high_lat)[0] # obtaining the parameters
    print consts_high_lat
    
    line_low = consts_low_lat[0]*np.arange(1.,9.)**(1./5.) + consts_low_lat[1] # regression line
    line_high = consts_high_lat[0]*np.arange(10.,31.)**(3./4.) + consts_high_lat[1] # regression line
    plt.plot(np.arange(1.,9.), line_low, 'r-')
    plt.plot(np.arange(10.,31.), line_high, 'b-')
    
    plt.plot(edge_ll_summer,psi_min_summer,'kx')
    ax = plt.gca()
    plt.xlabel('Northernmost latitude of SH Hadley cell')
    plt.ylabel('Peak strength of SH Hadley circulation')
    plt.xscale('log', subsx=[2,3,4,5])
    plt.yscale('log', subsx=[2,3,4,5])
    plt.xlim(1,30)
    plt.ylim(50,500)
    plt.xticks([1,5,10,20,30])
    ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
    plt.tight_layout()  
    figname = 'psi_latvsstrength_ylog.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    #plt.show()
    
    return edge_ll_summer, psi_min_summer


#psi_plot('full_qflux',lonin=[60.,150.])
#edge_ll_ap, psi_min_ap = psi_plot('ap_20')
edge_ll_ap, psi_min_ap = psi_plot('ap_2')
#edge_ll_fq, psi_min_fq = psi_plot('full_qflux')

#plt.plot(edge_ll_ap,psi_min_ap,'kx')
#ax = plt.gca()
#plt.plot(edge_ll_fq,psi_min_fq,'rx')
#plt.xlabel('Northernmost latitude of SH Hadley cell')
#plt.ylabel('Peak strength of SH Hadley circulation')
#plt.xscale('log', subsx=[2,3,4,5])
#plt.xlim(1,30)
#plt.xticks(xtick_vals)
#ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
#plt.show()
