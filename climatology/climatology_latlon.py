# Produce plots of the zonal mean terms in the momentum budget as a function of time for the aquaplanet simulations, showing the transition and change in balance of terms

from data_handling import time_means, season_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

    
def plot_clim(inp_fol):
    
    plot_dir = '/scratch/rg419/plots/climatology/'+inp_fol+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    sn_dic = season_dic(1)
    
    data = time_means(inp_fol,[121,181],filename='atmos_pentad',timeav='season')
    rain = (data.convection_rain + data.condensation_rain)*86400.
    del data
    
    data_p = time_means(inp_fol,[121,181],filename='atmos_pentad',timeav='season')
    ucomp = data_p.ucomp
    vcomp = data_p.vcomp
    del data_p
    
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/summer_holiday/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', ucomp.lat), ('lon', ucomp.lon)])
    
    for i in range(0,4):
        print i
        rain[i,:,:].plot.contourf(x='lon', y='lat', levels = np.arange(0.,31.,3.), add_label = False)
    #plt.ylim(-30,60)
    #plt.xlim(60,180)
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False)
        plt.quiver(ucomp.lon[::3], ucomp.lat[::1], ucomp[i,36,::1,::3], vcomp[i,36,::1,::3], scale=500.,headwidth=5)
        plt.savefig(plot_dir+'wind_and_rain_'+ sn_dic[i] + '.png')
        plt.close()    


plot_clim('era_amipsst')




