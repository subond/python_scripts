from data_handling import load_year_xr
from pentads import pentad_dic, month_dic
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

def load_phi(inp_fol,years,tper,trange=[0,0]):
    #set run name
    run_fol = inp_fol + '/np' + nproc
    year = years[0]
    rundata = load_year_xr(run_fol, year, pinterp=True)
    if tper == 'month':
        tvals = [30,11]
        tsize = 12
    elif tper == 'pentad':
        tvals = [5,71]
        tsize = trange[1]-trange[0]
    else:
        print 'time period not set'

    rundata.coords['tper'] = (rundata.time // tvals[0]) - tvals[1]
    mngrp = rundata.ucomp.groupby('tper').mean(('time'))
    if tper == 'pentad':
        dim = mngrp.tper.values[trange[0]:trange[1]]
    else:
        dim = mngrp.tper
    
    #Initialise arrays to load into
    phi =   xr.DataArray(np.zeros((tsize,17,64,128,len(years))), [(tper, dim ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    for year in years:
        print year
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords[tper] = (rundata.time // tvals[0]) - tvals[1]
        if tper == 'pentad':
            phi[:,:,:,:,year-years[0]]   = rundata.height.groupby(tper).mean(('time'))[trange[0]:trange[1],:,:,:]    
        else:
            phi[:,:,:,:,year-years[0]]   = rundata.height.groupby(tper).mean(('time'))    
    return phi.mean(('year'))/1000.


def plot_pd_gp_fn(inp_fol, nproc, years, trange):
    
    phi=load_phi(inp_fol, years, 'pentad', trange)
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/flat_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', phi.lat), ('lon', phi.lon)])
    pd_dic = pentad_dic(1)
    
    def plot_gp_var(lev,levels):
        g=phi[:,lev,:,:].plot.contourf(x='lon', y='lat',levels=levels, col='pentad', col_wrap=5)
        plt.ylim(0,45)
        plt.xlim(60,180)
        for i, ax in enumerate(g.axes.flat):
            if not 'aquaplanet' in inp_fol:
                land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
            ax.set_title(pd_dic[trange[0]+i+1])
        plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd_phi_'+ str(int(phi.pfull[lev].values)) +'.png')
        plt.close()
        
    plot_gp_var(2,np.arange(1.,1.6,0.025))
    plot_gp_var(9,np.arange(12.,14.1,0.1))            


def plot_mn_gp_fn(inp_fol, nproc, years):
    
    phi=load_phi(inp_fol, years, 'month')
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/flat_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', phi.lat), ('lon', phi.lon)])
    mn_dic = month_dic(1)

    def plot_gp_var(lev,levels):
        g=phi[:,lev,:,:].plot.contourf(x='lon', y='lat',levels=levels, col='month', col_wrap=4)
        plt.ylim(0,45)
        plt.xlim(60,180)
        for i, ax in enumerate(g.axes.flat):
            if not 'aquaplanet' in inp_fol:
                land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
            ax.set_title(mn_dic[i+1])
        plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_phi_'+ str(int(phi.pfull[lev].values)) +'.png')
        plt.close()
        
    plot_gp_var(2,np.arange(0.,2.1,0.1))
    plot_gp_var(9,np.arange(12.,14.1,0.1))
        

    
if __name__ == "__main__":

    nproc   = '16'

    #plot_mn_gp_fn('aquaplanet_1m_diags', nproc,range(2,17))
    #plot_mn_gp_fn('aquaplanet_10m_diags', nproc,range(2,17))
    #plot_mn_gp_fn('topo_10m_diags', nproc, range(2,17))
    #plot_mn_gp_fn('flat_10m_diags', nproc, range(2,17))
    
    #plot_pd_gp_fn('aquaplanet_1m_diags', nproc,range(2,17),[25,45])
    #plot_pd_gp_fn('aquaplanet_10m_diags', nproc,range(2,17),[40,60])
    plot_pd_gp_fn('topo_10m_diags', nproc, range(2,17),[35,55])
    plot_pd_gp_fn('flat_10m_diags', nproc, range(2,17),[35,55])
    

