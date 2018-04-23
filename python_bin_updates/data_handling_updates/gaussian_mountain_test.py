# Test how varying parameters gives different Gaussian mountain shapes

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os


# specify resolution
t_res = 42
#read in grid from approriate file
GFDL_BASE = os.environ['GFDL_BASE']
resolution_file = Dataset(GFDL_BASE + 'src/extra/python/scripts/gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')
lons = resolution_file.variables['lon'][:]
lats = resolution_file.variables['lat'][:]
lonb = resolution_file.variables['lonb'][:]
latb = resolution_file.variables['latb'][:]
nlon=lons.shape[0]
nlat=lats.shape[0]
#make 2d arrays of latitude and longitude
lon_array, lat_array = np.meshgrid(lons,lats)
lonb_array, latb_array = np.meshgrid(lonb,latb)
    
    
def mountain_test(h_0=2670., central_lon=247.5, central_lat=40, l_1=7.5, l_2=20., gamma_1=42., gamma_2=42.):
    
    topo_array = np.zeros((nlat,nlon))
    
    delta_1 = ((lon_array - central_lon)*np.cos(np.radians(gamma_1)) + (lat_array - central_lat)*np.sin(np.radians(gamma_1)))/l_1
    delta_2 = (-(lon_array - central_lon)*np.sin(np.radians(gamma_2)) + (lat_array - central_lat)*np.cos(np.radians(gamma_2)))/l_2
    h_arr = h_0 * np.exp(-(delta_1**2. + delta_2**2.))
    idx = (h_arr / h_0 > 0.05) 
    
    topo_array[idx] = h_arr[idx]
    
    return topo_array


def plot_check(topo_array, fig_no=1):
    #Show configuration on screen to allow checking
    lon_0 = lons.mean()
    lat_0 = lats.mean()
    m = Basemap(lat_0=lat_0,lon_0=lon_0)
    xi, yi = m(lon_array, lat_array)
    plt.figure(fig_no)
    cs = m.contourf(xi,yi,topo_array, cmap=plt.get_cmap('RdBu_r'))
    cb = plt.colorbar(cs, shrink=0.5, extend='both')
    plt.xticks(np.linspace(0,360,13))
    plt.yticks(np.linspace(-90,90,7))



if __name__ == "__main__":

    topo_array_1 = mountain_test()
    plot_check(topo_array_1, fig_no=1)
    
    topo_array_2 = mountain_test(gamma_1=8., gamma_2=8., central_lon=285., central_lat=-19., l_1=3.75, l_2=36, h_0=4080.)
    plot_check(topo_array_2, fig_no=2)
    
    topo_array_3 = mountain_test(gamma_1=-36., gamma_2=-36., central_lon=43., central_lat=-6., l_1=7.5, l_2=36, h_0=1530.)
    plot_check(topo_array_3, fig_no=3)
    
    
    
    plt.show()

