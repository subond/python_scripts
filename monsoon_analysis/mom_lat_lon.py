# Make contour plots of the terms in the momentum budget at 200hPa plus the wind vectors for the topography run over the onset period. 

from physics import mombudg_lev_pd_fn
from pentads import pentad_dic, month_dic
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


def plot_pd_mom_fn(inp_fol, years, trange):
    
    data=mombudg_lev_pd_fn(inp_fol, years, trange)
    print 'data loaded, starting plotting'
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/topo_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
    pd_dic = pentad_dic(1)
    
    data.to_netcdf(path='/scratch/rg419/plots/monsoon_analysis/'+inp_fol+'/mombudg_data.nc', mode='w')
    
    def plot_mom_var(var,levels,qscale):
        
        if var == 'fv_mn_imb':
            plot_data = data.data_vars['fv_av'] + data.data_vars['mom_mean']
        else:
            plot_data = data.data_vars[var]
        plot_data = plot_data*10000.
        
        for mnth in range(0,12):
            g=plot_data[mnth*6 : mnth*6+6,:,:].plot.contourf(x='lon', y='lat',levels=levels, col='pentad', col_wrap=3)
            plt.ylim(-45,45)
            plt.xlim(60,180)
            for i, ax in enumerate(g.axes.flat):
                land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
                ax.quiver(data.lon[::3], data.lat[::1], data.u_av[mnth*6 +i,::1,::3], data.v_av[mnth*6+i,::1,::3], scale=qscale,headwidth=5)
                ax.set_title(pd_dic[mnth*6+i+1])
            plt.savefig('/scratch/rg419/plots/monsoon_analysis/'+ inp_fol + '/' + var + str(mnth + 1) + '.png')
            plt.close()
        
    vars = ['dphidx_av','fv_av','mom_eddy','mom_mean','mom_sum','fv_mn_imb']
    for var in vars:
        plot_mom_var(var,np.arange(-2.,2.1,0.2),500.)            


plot_pd_mom_fn('aquaplanet_10m',range(11,41),[0,72])
print 'aquaplanet_10m done'
#plot_pd_mom_fn('aquaplanet_2m',range(11,41),[0,72])
#print 'aquaplanet_2m done'
#plot_pd_mom_fn('flat_10m',range(11,41),[0,72])
#print 'flat_10m done'
#plot_pd_mom_fn('aquamountain_10m',range(11,41),[0,72])
#print 'aquamountain_10m done'
#plot_pd_mom_fn('topo_10m',range(11,41),[0,72])
#print 'topo_10m done'
#plot_pd_mom_fn('flat_10m',range(21,41),[31,51])
