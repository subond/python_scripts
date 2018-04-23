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

#import statsmodels.api as sm


    
rcParams['figure.figsize'] = 9, 6
rcParams['font.size'] = 20
rcParams['text.usetex'] = True
    
plot_dir = '/scratch/rg419/plots/rotation/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

lev = 500.
thresh = 0.

def get_edge_psi(run, lonin=[-1.,361.], sanity_check=False, lev=lev):
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(data.lat * np.pi/180.)
    
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


def regime_plot(run, ylims=[100,600], yticks=[100,200,300,400,500,600]):
    edge_ll, psi_min = get_edge_psi(run)
    
    plt.plot(edge_ll, psi_min, 'kx', ms=10, mew=2)
    plt.plot(edge_ll[40:50], psi_min[40:50], 'rx', ms=10, mew=2)
    ax = plt.gca()
    plt.xlabel('$\phi_{0}$')
    plt.ylabel('$\Psi_{max}$, 10$^9$ kg/s')
    plt.xscale('log', subsx=[2,3,4,5])
    plt.yscale('log', subsx=[2,3,4,5])
    plt.xlim(1,40)
    plt.ylim(ylims[0],ylims[1])
    plt.xticks([1,5,10,20,40])
    plt.yticks(yticks)
    ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(tk.ScalarFormatter())
    plt.tight_layout()  
    
    figname = 'psi_latvsstrength_' + run + '_' + str(int(lev)) + '_' + str(int(thresh)) + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()

regime_plot('sn_1.000')
regime_plot('rt_0.500', ylims=[200,600], yticks=[200,300,400,500,600])
regime_plot('rt_2.000', ylims=[40,240], yticks=[40,80,120,160,200,240])
