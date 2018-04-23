from momentum_budget_fns import mombudg_2d_pd_fn, mombudg_2d_mn_fn, mombudg_2d_an_fn
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


def monsoon_onset_mom_fn(inp_fol, nproc, range):
    
    data=mombudg_2d_pd_fn(inp_fol, nproc, years)
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/flat_10m_diags/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
    
    g = data.dphidx_av[range[0]:range[1],9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_dphidx_200.png')
    plt.close()
    g = data.fv_av[range[0]:range[1],9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_fv_200.png')
    plt.close()
    g = data.mom_eddy[range[0]:range[1],9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_momeddy_200.png')
    plt.close()
    g = data.mom_mean[range[0]:range[1],9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_mommean_200.png')
    plt.close()
    g = data.ddamp_av[range[0]:range[1],9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_ddamp_200.png')
    plt.close()
    g = data.mom_sum[range[0]:range[1],9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_momsum_200.png')
    plt.close()

    g = data.dphidx_av[range[0]:range[1],2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_dphidx_850.png')
    plt.close()
    g = data.fv_av[range[0]:range[1],2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_fv_850.png')
    plt.close()
    g = data.mom_eddy[range[0]:range[1],2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_momeddy_850.png')
    plt.close()
    g = data.mom_mean[range[0]:range[1],2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_mommean_850.png')
    plt.close()
    g = data.ddamp_av[range[0]:range[1],2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_ddamp_850.png')
    plt.close()
    g = data.mom_sum[range[0]:range[1],2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='pentad', col_wrap=5)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/pd2_momsum_850.png')
    plt.close()



def plot_mn_mom_fn(inp_fol, nproc, years):
    
    data=mombudg_2d_mn_fn(inp_fol, nproc, years)
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/flat_10m_diags/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])

    g = data.dphidx_av[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_dphidx_200.png')
    plt.close()
    g = data.fv_av[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_fv_200.png')
    plt.close()
    g = data.mom_eddy[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_momeddy_200.png')
    plt.close()
    g = data.mom_mean[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_mommean_200.png')
    plt.close()
    g = data.ddamp_av[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_ddamp_200.png')
    plt.close()
    g = data.mom_sum[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_momsum_200.png')
    plt.close()

    g = data.dphidx_av[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_dphidx_850.png')
    plt.close()
    g = data.fv_av[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_fv_850.png')
    plt.close()
    g = data.mom_eddy[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_momeddy_850.png')
    plt.close()
    g = data.mom_mean[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_mommean_850.png')
    plt.close()
    g = data.ddamp_av[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_ddamp_850.png')
    plt.close()
    g = data.mom_sum[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001), col='month', col_wrap=4)
    plt.ylim(0,45)
    plt.xlim(60,180)
    for i, ax in enumerate(g.axes.flat):
        land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/mn_momsum_850.png')
    plt.close()



def plot_an_mom_fn(inp_fol, nproc, years):
    
    data=mombudg_2d_an_fn(inp_fol, nproc, years)
    
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/flat_10m_diags/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])

    data.dphidx_av[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_dphidx_200.png')
    plt.close()
    data.fv_av[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_fv_200.png')
    plt.close()
    data.mom_eddy[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_momeddy_200.png')
    plt.close()
    data.mom_mean[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_mommean_200.png')
    plt.close()
    data.ddamp_av[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_ddamp_200.png')
    plt.close()
    data.mom_sum[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_momsum_200.png')
    plt.close()

    data.dphidx_av[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_dphidx_850.png')
    plt.close()
    data.fv_av[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_fv_850.png')
    plt.close()
    data.mom_eddy[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_momeddy_850.png')
    plt.close()
    data.mom_mean[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_mommean_850.png')
    plt.close()
    data.ddamp_av[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_ddamp_850.png')
    plt.close()
    data.mom_sum[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False)
    #plt.ylim(0,45)
    #plt.xlim(60,180)
    plt.savefig('/scratch/rg419/plots/'+inp_fol+'/an_momsum_850.png')
    plt.close()


    
if __name__ == "__main__":

    nproc   = '16'

    plot_mn_mom_fn('flat_10m_diags', nproc, range(2,17))
    #plot_mn_mom_fn('topo_10m', nproc, range(2,7))

    #monsoon_onset_mom_fn('aquaplanet_1m_diags', nproc, [25,45])
    #print '1m done'
    #monsoon_onset_mom_fn('aquaplanet_10m_diags', nproc, [40,60])
    #print '10m done'
    #monsoon_onset_mom_fn('flat_10m_diags', nproc, [35,55])
    #print 'All done :)'
