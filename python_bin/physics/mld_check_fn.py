#Load up mixed layer depth and check values are correct

from netCDF4 import Dataset
import sys
sys.path.append('/scratch/rg419/workdir_moist/python_bin')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def mld_check_fn(run_fol, inp_fol):

    #load up gridding details
    nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/' + run_fol + '/run1/atmos_daily.nc'
    fh = Dataset(nc_file, mode='r')
    lats = fh.variables['lat'][:]
    lons = fh.variables['lon'][:]
    mld = fh.variables['ml_heat_cap'][:]
    fh.close()

    rho_cp = 1.035e3 * 3989.24495292815

    plt.set_cmap('jet')
    plt.contourf(lons,lats,mld/rho_cp, np.arange(0,30,1))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.colorbar()
    plt.savefig('/scratch/rg419/workdir_moist/' + inp_fol + '/mld.png')
    plt.clf()

    print np.amin(np.amin(mld/rho_cp)), np.amax(np.amax(mld/rho_cp))

    return mld