#Produce plots comparable to Bordoni and Schneider 2008 to test mixed layer depth effects
#need precip averaged over Asia - use 70 to 130. Compare also with other regions.

from netCDF4 import Dataset
import sys
sys.path.append('/scratch/rg419/workdir_moist/python_bin')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from plotting_fns import plot_var_ll, plot_quiver_ll, plot_var_zav
from time_av_fns import multi_month_load, monthly_mean, intermonthly_mean

def bs_plot_fn(run_fol, inp_fol, years):

    #load up gridding details
    nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/' + run_fol + '/run1/atmos_daily.nc'
    fh = Dataset(nc_file, mode='r')
    lats = fh.variables['lat'][:]
    lons = fh.variables['lon'][:]
    fh.close()

    def load_year_precip(year):
        #read in a year's worth of precip data
        cnvp = multi_month_load('convection_rain', run_fol, range(year*12+1,year*12+13), twodim=True)
        cndp = multi_month_load('condensation_rain', run_fol, range(year*12+1,year*12+13), twodim=True)
        totp = (cnvp+cndp)*86400

        #take zonal mean over asia
        #print lons[25], lons[46]   #display latitude range
        totp_asia = np.mean(totp[:,:,:,25:46],3)
        #reshape time dimensions from month,day to day
        totp_asia_ts = np.rollaxis( np.reshape( totp_asia ,(360,64)),1)
        #take pentad average
        totp_asia_pd = np.zeros((64,72))
        for i in range(0,72):
            totp_asia_pd[:,i] = np.mean(totp_asia_ts[:,i*5:i*5+4],1)

        #take global mean to compare
        totp_glob = np.mean(totp,3)
        totp_glob_ts = np.rollaxis( np.reshape( totp_glob ,(360,64)),1)
        totp_glob_pd = np.zeros((64,72))
        for i in range(0,72):
            totp_glob_pd[:,i] = np.mean(totp_glob_ts[:,i*5:i*5+4],1)

        return totp_asia_pd, totp_glob_pd
    
    #load in data for requested years
    totp_asia_yrs = np.zeros((64,72,len(years)))
    totp_glob_yrs = np.zeros((64,72,len(years)))
    for year in years:
        totp_asia_yrs[:,:,year-1], totp_glob_yrs[:,:,year-1] = load_year_precip(year)

    #average over years
    totp_asia_pd = np.mean(totp_asia_yrs,2)
    totp_glob_pd = np.mean(totp_glob_yrs,2)

    #plot and save 
    pd=range(3,361,5)
    plt.set_cmap('hot_r')
    plt.contourf(pd,lats,totp_asia_pd,range(2,21,2))
    plt.ylim((-40,40))
    plt.xlabel('Day')
    plt.ylabel('Latitude')
    plt.colorbar()
    plt.savefig('/scratch/rg419/workdir_moist/' + inp_fol + '/totp_asia.png')
    plt.clf()

    plt.set_cmap('hot_r')
    plt.contourf(pd,lats,totp_glob_pd,range(2,21,2))
    plt.ylim((-40,40))
    plt.xlabel('Day')
    plt.ylabel('Latitude')
    plt.colorbar()
    plt.savefig('/scratch/rg419/workdir_moist/' + inp_fol + '/totp_glob.png')
    plt.clf()

    return totp_asia_pd, totp_glob_pd




def itcz_lat_fn(run_fol, inp_fol, years):

    #read in lats
    nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/' + run_fol + '/run1/atmos_daily.nc'
    fh = Dataset(nc_file, mode='r')
    lats = fh.variables['lat'][:]
    fh.close()

    def find_itcz_lat(year):
        #load up surface temperatures as latitude time array
        t_surf_in = multi_month_load('t_surf', run_fol, range(year*12+1,year*12+13), twodim=True)
        tsurf = np.reshape(np.mean(t_surf_in,3),(360,64))

        #locate peak in zonal mean surface temperature
        itcz_lat = np.zeros(360)
        for i in range(0,360):
            tmax = np.amax(tsurf[i,:])
            itcz_lat[i] = lats[[k for k, j in enumerate(tsurf[i,:]) if j == tmax]]
        return itcz_lat

    #load in data for requested years
    itcz_lat_yrs = np.zeros((360,len(years)))
    for year in years:
        itcz_lat_yrs[:,year-1] = find_itcz_lat(year)

    #average over years
    itcz_lat = np.mean(itcz_lat_yrs,1)

    #plot and save scatter
    plt.plot(itcz_lat,'x',markersize=3)
    plt.ylim((-30,30))
    plt.xlabel('Day')
    plt.ylabel('Peak SST latitude')
    plt.savefig('/scratch/rg419/workdir_moist/' + inp_fol + '/itcz_lat.png')
    plt.clf()

    return itcz_lat
