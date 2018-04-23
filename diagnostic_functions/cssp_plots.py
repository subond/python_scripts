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

def plot_mn(mrange,dir):

    m_dic = {13:'jan', 14:'feb', 15:'mar', 16:'apr', 17:'may', 18:'june', 19:'jul', 20:'aug', 21:'sept', 22:'oct', 23:'nov', 24:'dec' }

    u = monthly_mean('ucomp', run_fol, mrange, twodim=False)
    v = monthly_mean('vcomp', run_fol, mrange, twodim=False)
    isotach = np.sqrt(u**2 + v**2)
    w = monthly_mean('omega', run_fol, mrange, twodim=False)
    temp = monthly_mean('temp', run_fol, mrange, twodim=False)
    tsurf = monthly_mean('t_surf', run_fol, mrange, twodim=True)

    flux_t = -monthly_mean('flux_t', run_fol, mrange, twodim=True)
    flux_sw = monthly_mean('flux_sw', run_fol, mrange, twodim=True)
    flux_lhe = -monthly_mean('flux_lhe', run_fol, mrange, twodim=True)
    flux_lw = monthly_mean('flux_lw', run_fol, mrange, twodim=True) - stefan*tsurf**4

    cnvp = monthly_mean('convection_rain',run_fol, mrange,twodim=True)
    cndp = monthly_mean('condensation_rain',run_fol, mrange,twodim=True)
    totp = (cnvp+cndp)*86400

    for m in mrange:
        plt.set_cmap('bwr')
        plot_var_ll(lons,lats,tsurf[m-13,:,:], v=range(230,321,5), plotno=1,land_mask=land_mask, pltitle = 'Surface temperature ('+m_dic[m]+')', cblabel = 'K')
        plt.savefig(dir + 'tsurf_'+str(m)+'.png')
        plt.clf()

        #plt.set_cmap('bwr')
        plot_var_ll(lons,lats,w[m-13,37,:,:], v=np.arange(-0.2,0.21,0.01), plotno=1, land_mask=land_mask, pltitle = 'Vertical velocity at 710 hPa ('+m_dic[m]+')', cblabel = 'Pa/s')
        plt.savefig(dir + 'w710_'+str(m)+'.png')
        plt.clf()

        plt.set_cmap('viridis_r')
        plot_var_ll(lons,lats,isotach[m-13,38,:,:], v=np.arange(2,32,1), plotno=1, land_mask=land_mask, cblabel = 'm/s')
        plot_quiver_ll(lons[::2],lats[::2],u[m-13,38,::2,::2], v[m-13,38,::2,::2], plotno=1,  pltitle = 'Vector wind with isotachs at 815 hPa ('+m_dic[m]+')' ,draw_grid = False)
        plt.savefig(dir + 'uv815_'+str(m)+'.png')
        plt.clf()

        plt.set_cmap('RdYlBu')
        plot_var_ll(lons,lats,flux_t[m-13,:,:], v=range(-60,61,5), plotno=1,land_mask=land_mask, pltitle = 'Sensible heat flux, positive downwards ('+m_dic[m]+')', cblabel = 'W/m2')
        plt.savefig(dir + 'flux_t'+str(m)+'.png')
        plt.clf()

        plt.set_cmap('autumn')
        plot_var_ll(lons,lats,flux_lhe[m-13,:,:], plotno=1, v=range(-280,0,10), land_mask=land_mask, pltitle = 'Latent heat flux, positive downwards ('+m_dic[m]+')', cblabel = 'W/m2')
        plt.savefig(dir + 'flux_lhe'+str(m)+'.png')
        plt.clf()

        plt.set_cmap('RdYlBu_r')
        plot_var_ll(lons,lats,flux_sw[m-13,:,:], v=range(0,321,10), plotno=1,land_mask=land_mask, pltitle = 'Surface solar radiation flux, positive downwards ('+m_dic[m]+')', cblabel = 'W/m2')
        plt.savefig(dir + 'flux_sw'+str(m)+'.png')
        plt.clf()

        plt.set_cmap('viridis')
        plot_var_ll(lons,lats,flux_lw[m-13,:,:] , v=range(-130,-39,5), plotno=1,land_mask=land_mask, pltitle = 'Surface thermal radiation flux, positive downwards ('+m_dic[m]+')', cblabel = 'W/m2')
        plt.savefig(dir + 'flux_lw'+str(m)+'.png')
        plt.clf()

        plt.set_cmap('Blues')
        plot_var_ll(lons,lats,totp[m-13,:,:], v=range(0,16,1), plotno=1,land_mask=land_mask, pltitle = 'Total precipitation ('+m_dic[m]+')', cblabel = 'mm/day')
        plt.savefig(dir + 'totp'+str(m)+'.png')
        plt.clf()


    return

plot_mn(range(13,25),'monsoonplot/')



