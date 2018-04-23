from physics.momentum_budget_fns import mombudg_2d_pd_fn, mombudg_2d_mn_fn, mombudg_2d_an_fn
from pentads import pentad_dic, month_dic
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


def plot_pd_mom_fn(inp_fol, nproc, years, trange):
    
    data=mombudg_2d_pd_fn(inp_fol, nproc, years, trange)
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/flat_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
    pd_dic = pentad_dic(1)
    
    def plot_mom_var(var,lev,levels,qscale):
        g=data.data_vars[var][:,lev,:,:].plot.contourf(x='lon', y='lat',levels=levels, col='pentad', col_wrap=5)
        plt.ylim(0,45)
        plt.xlim(60,180)
        for i, ax in enumerate(g.axes.flat):
            if not 'aquaplanet' in inp_fol:
                land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
            ax.quiver(data.lon[::3], data.lat[::1], data.u_av[i,lev,::1,::3], data.v_av[i,lev,::1,::3], scale=qscale,headwidth=5)
            ax.set_title(pd_dic[trange[0]+i+1])
        plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd_'+var+'_'+ str(int(data.pfull[lev].values)) +'.png')
        plt.close()
        
    vars = ['dphidx_av','fv_av','mom_eddy','mom_mean','ddamp_av','mom_sum']
    for var in vars:
        plot_mom_var(var,2,np.arange(-0.0003,0.00031,0.00003),250.)
        plot_mom_var(var,9,np.arange(-0.0005,0.00051,0.00005),500.)            

def plot_mn_mom_fn(inp_fol, nproc, years):
    
    data=mombudg_2d_mn_fn(inp_fol, nproc, years)
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/flat_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
    mn_dic = month_dic(1)

    def plot_mom_var(var,lev,levels,qscale):
        g=data.data_vars[var][:,lev,:,:].plot.contourf(x='lon', y='lat',levels=levels, col='month', col_wrap=4)
        plt.ylim(0,45)
        plt.xlim(60,180)
        for i, ax in enumerate(g.axes.flat):
            if not 'aquaplanet' in inp_fol:
                land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
            ax.quiver(data.lon[::3], data.lat[::1], data.u_av[i,lev,::1,::3], data.v_av[i,lev,::1,::3], scale=qscale,headwidth=5)
            ax.set_title(mn_dic[i+1])
        plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_'+var+'_'+ str(int(data.pfull[lev].values)) +'.png')
        plt.close()
        
    vars = ['dphidx_av','fv_av','mom_eddy','mom_mean','ddamp_av','mom_sum']
    for var in vars:
        plot_mom_var(var,2,np.arange(-0.0003,0.00031,0.00003),250.)
        plot_mom_var(var,9,np.arange(-0.0005,0.00051,0.00005),500.)
        

def plot_an_mom_fn(inp_fol, nproc, years):
    
    data=mombudg_2d_an_fn(inp_fol, nproc, years)  
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/flat_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])

    def plot_mom_var(var,lev,levels):
        data.data_vars[var][:,lev,:,:].plot.contourf(x='lon', y='lat',levels=levels)
        if not 'aquaplanet' in inp_fol:
            land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False)
        plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_'+var+'_'+ str(int(data.pfull[lev].values)) +'.png')
        plt.close()
        
    vars = ['dphidx_av','fv_av','mom_eddy','mom_mean','ddamp_av','mom_sum']
    for var in vars:
        for lev in range(2,10,7):
            plot_mom_var(var,lev,np.arange(-0.0002,0.0002,0.00001))


    
if __name__ == "__main__":

    nproc   = '16'

    #plot_mn_mom_fn('aquaplanet_1m_diags', nproc,range(2,17))
    #plot_mn_mom_fn('aquaplanet_10m_diags', nproc,range(2,17))
    #plot_mn_mom_fn('topo_10m', nproc, range(2,17))
    #plot_mn_mom_fn('flat_10m_diags', nproc, range(2,17))
    
    #plot_pd_mom_fn('aquaplanet_1m_diags', nproc,range(2,17),[25,45])
    #plot_pd_mom_fn('aquaplanet_10m_diags', nproc,range(2,17),[40,60])
    #plot_pd_mom_fn('topo_10m', nproc, range(2,17),[35,55])
    plot_pd_mom_fn('flat_10m_diags', nproc, range(2,17),[35,55])
    

