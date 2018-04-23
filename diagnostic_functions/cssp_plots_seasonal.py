from netCDF4 import Dataset

import sys
sys.path.append('/scratch/rg419/workdir_moist')

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from plotting_fns import plot_var_ll, plot_quiver_ll, plot_var_zav
from time_av_fns import multi_month_load, monthly_mean, intermonthly_mean

#set run name
run_fol = 'cssp_run_topo/np16'
inp_fol = 'cssp_run_topo'

#run_fol = 'land_test/np16'
#inp_fol = 'land_test'

#open files

nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/' + run_fol + '/run1/atmos_daily.nc'
fh = Dataset(nc_file, mode='r')
land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/' + inp_fol + '/input/land.nc'
lfh = Dataset(land_file, 'r', format='NETCDF3_CLASSIC')
land_mask = lfh.variables['land_mask'][:]
lfh.close()

#read in variables
lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
p = fh.variables['pfull'][:]
times = fh.variables['time'][:]
fh.close()

stefan = 5.6734e-8

#lons = np.roll(lons,-22)
#print lons

def plot_season(season,dir):

    year_dic = {'djf': [13,14,24], 'mam': [15,16,17], 'jja': [18,19,20], 'son': [21,22,23]}

    u = intermonthly_mean('ucomp', run_fol, year_dic[season], twodim=False)
    v = intermonthly_mean('vcomp', run_fol, year_dic[season], twodim=False)
    isotach = np.sqrt(u**2 + v**2)
    w = intermonthly_mean('omega', run_fol, year_dic[season], twodim=False)
    temp = intermonthly_mean('temp', run_fol, year_dic[season], twodim=False)
    tsurf = intermonthly_mean('t_surf', run_fol, year_dic[season], twodim=True)

    flux_t = -intermonthly_mean('flux_t', run_fol, year_dic[season], twodim=True)
    flux_sw = intermonthly_mean('flux_sw', run_fol, year_dic[season], twodim=True)
    flux_lhe = -intermonthly_mean('flux_lhe', run_fol, year_dic[season], twodim=True)
    flux_lw = intermonthly_mean('flux_lw', run_fol, year_dic[season], twodim=True) - stefan*tsurf**4

    cnvp = intermonthly_mean('convection_rain',run_fol, year_dic[season],twodim=True)
    cndp = intermonthly_mean('condensation_rain',run_fol, year_dic[season],twodim=True)
    totp = (cnvp+cndp)*86400


    plt.set_cmap('bwr')
    plot_var_ll(lons,lats,tsurf[:,:], v=range(230,321,5), plotno=1,land_mask=land_mask, pltitle = 'Surface temperature ('+season+')', cblabel = 'K',lon_0=[60.,240.], cbarloc='right')
    plt.savefig(dir + 'tsurf_'+season+'.png')
    plt.clf()

    #plt.set_cmap('bwr')
    plot_var_ll(lons,lats,w[37,:,:], v=np.arange(-0.2,0.21,0.01), plotno=1, land_mask=land_mask, pltitle = 'Vertical velocity at 710 hPa ('+season+')', cblabel = 'Pa/s',lon_0=[60.,240.], cbarloc='right')
    plt.savefig(dir + 'w710_'+season+'.png')
    plt.clf()

    plt.set_cmap('viridis_r')
    plot_var_ll(lons,lats,isotach[38,:,:], v=np.arange(2,32,1), plotno=1, land_mask=land_mask, cblabel = 'm/s',lon_0=[60.,240.], cbarloc='right')
    plot_quiver_ll(lons[::2],lats[::2],u[38,::2,::2], v[38,::2,::2], plotno=1,  pltitle = 'Vector wind with isotachs at 815 hPa ('+season+')' ,draw_grid = False,lon_0=[60.,240.])
    plt.savefig(dir + 'uv815_'+season+'.png')
    plt.clf()

    plt.set_cmap('RdYlBu')
    plot_var_ll(lons,lats,flux_t[:,:], v=range(-60,61,5), plotno=1,land_mask=land_mask, pltitle = 'Sensible heat flux, positive downwards ('+season+')', cblabel = 'W/m2',lon_0=[60.,240.], cbarloc='right')
    plt.savefig(dir + 'flux_t_'+season+'.png')
    plt.clf()

    plt.set_cmap('autumn')
    plot_var_ll(lons,lats,flux_lhe[:,:], plotno=1, v=range(-280,0,10), land_mask=land_mask, pltitle = 'Latent heat flux, positive downwards ('+season+')', cblabel = 'W/m2',lon_0=[60.,240.], cbarloc='right')
    plt.savefig(dir + 'flux_lhe_'+season+'.png')
    plt.clf()

    plt.set_cmap('RdYlBu_r')
    plot_var_ll(lons,lats,flux_sw[:,:], v=range(0,321,10), plotno=1,land_mask=land_mask, pltitle = 'Surface solar radiation flux, positive downwards ('+season+')', cblabel = 'W/m2',lon_0=[60.,240.], cbarloc='right')
    plt.savefig(dir + 'flux_sw_'+season+'.png')
    plt.clf()

    plt.set_cmap('viridis')
    plot_var_ll(lons,lats,flux_lw[:,:] , v=range(-130,-39,5), plotno=1,land_mask=land_mask, pltitle = 'Surface thermal radiation flux, positive downwards ('+season+')', cblabel = 'W/m2',lon_0=[60.,240.], cbarloc='right')
    plt.savefig(dir + 'flux_lw_'+season+'.png')
    plt.clf()

    plt.set_cmap('Blues')
    plot_var_ll(lons,lats,totp[:,:], v=range(0,16,1), plotno=1,land_mask=land_mask, pltitle = 'Total precipitation ('+season+')', cblabel = 'mm/day',lon_0=[60.,240.], cbarloc='right')
    plt.savefig(dir + 'totp_'+season+'.png')
    plt.clf()


    return


plot_season('djf','monsoonplot/')
plot_season('mam','monsoonplot/')
plot_season('jja','monsoonplot/')
plot_season('son','monsoonplot/')



