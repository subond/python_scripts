# Make contour plots of the terms in the momentum budget at 200hPa plus the wind vectors for the topography run over the onset period. Use 10 day averages, play with the choice of averaging period, and ideally have only 4 or so panels.

from physics import mombudg_lev_pd_fn
from pentads import pentad_dic, month_dic
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


def plot_pd_mom_fn(inp_fol, years, trange):
    
    data=mombudg_lev_pd_fn(inp_fol, years, trange)
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/topo_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
    pd_dic = pentad_dic(1)
    
    def plot_mom_var(var,levels,qscale):
        var_dic = {'fv_av': 'fv',
                    'mom_eddy': 'Eddy advective terms',
                    'mom_mean': 'Mean state advective terms',
                    'mom_sum': 'Residual',
                    'dphidx_av': 'Geopotential gradient'}
        plot_data = data.data_vars[var]
        plot_data = plot_data*10000.
        g=plot_data.plot.contourf(x='lon', y='lat',levels=levels, col='pentad', col_wrap=5)
        plt.ylim(-45,45)
        plt.xlim(60,180)
        for i, ax in enumerate(g.axes.flat):
            land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
            ax.quiver(data.lon[::3], data.lat[::1], data.u_av[i,::1,::3], data.v_av[i,::1,::3], scale=qscale,headwidth=5)
            ax.set_title(pd_dic[trange[0]+i+1])
        plt.savefig('/scratch/rg419/plots/mom_budg_work/'+ inp_fol + '/mom_ll_' + var + '.png')
        plt.close()
        
    vars = ['dphidx_av','fv_av','mom_eddy','mom_mean','mom_sum']
    for var in vars:
        plot_mom_var(var,np.arange(-5.,5.1,0.5),500.)            


def plot_ll_mom_fn(inp_fol, years, trange):
    
    data=mombudg_lev_pd_fn(inp_fol, years, trange)
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/topo_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
    pd_dic = pentad_dic(1)
    
    def plot_mom_var(var,levels,qscale):
        var_dic = {'fv_av': 'fv',
                    'mom_eddy': 'Eddy advective terms',
                    'mom_mean': 'Mean state advective terms',
                    'mom_sum': 'Residual',
                    'dphidx_av': 'Geopotential gradient'}
        plot_data = data.data_vars[var]
        plot_data = plot_data[0:4,:,:].mean(('pentad'))*10000.
        g=plot_data.plot.contourf(x='lon', y='lat',levels=levels)
        plt.ylim(-45,45)
        plt.xlim(60,180)
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False)
        plt.quiver(data.lon[::3], data.lat[::1], data.u_av[0:4,::1,::3].mean('pentad'), data.v_av[0:4,::1,::3].mean('pentad'), scale=qscale,headwidth=5)
        plt.savefig('/scratch/rg419/plots/mom_budg_work/'+ inp_fol + '/mom_ll_' + var + '_before.png')
        plt.close()
        
        plot_data = data.data_vars[var]
        plot_data = plot_data[5:9,:,:].mean(('pentad'))*10000.
        g=plot_data.plot.contourf(x='lon', y='lat',levels=levels)
        plt.ylim(-45,45)
        plt.xlim(60,180)
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False)
        plt.quiver(data.lon[::3], data.lat[::1], data.u_av[5:9,::1,::3].mean('pentad'), data.v_av[5:9,::1,::3].mean('pentad'), scale=qscale,headwidth=5)
        plt.savefig('/scratch/rg419/plots/mom_budg_work/'+ inp_fol + '/mom_ll_' + var + '_after.png')
        plt.close()
        
    vars = ['dphidx_av','fv_av','mom_eddy','mom_mean','mom_sum']
    for var in vars:
        plot_mom_var(var,np.arange(-5.,5.1,0.5),500.)            



plot_ll_mom_fn('topo_10m',range(21,41),[36,45])
#plot_pd_mom_fn('flat_10m',range(21,41),[31,51])
