# Produce plots of the zonal mean terms in the momentum budget as a function of time for the aquaplanet simulations, showing the transition and change in balance of terms

from data_handling import pentad_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#plt.rc('text', usetex=True)

def plot_pd_mom_fn(inp_fol):
    
    data_file = '/scratch/rg419/plots/monsoon_analysis/'+inp_fol+'/mom_stat_data.nc'
    data= xr.open_dataset( data_file)

    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/topo_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
    pd_dic = pentad_dic(1)
    
    #mom_sum = (data.dphidx*9.8 + data.fv - data.mom_trans - data.mom_stat - data.mom_mean)*10000.
    #mom_sum[60,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-2.,2.1,0.2))
    #plt.show()
        
    def plot_mom_var(var,levels,qscale):
        
        if var == 'mom_sum':
            plot_data = data.dphidx*9.8 + data.fv - data.mom_trans - data.mom_stat - data.mom_mean
        else:
            plot_data = data.data_vars[var]
        if var == 'dphidx':
            plot_data = plot_data*9.8
        if var == 'fv':
            plot_data = plot_data - data.mom_mean
        plot_data = plot_data*10000.
        #print plot_data.values.shape
        
        for mnth in range(0,12):
            if var == 'fv':
                g=plot_data[:,mnth*6 : mnth*6+6,9,:].plot.contourf(x='lon', y='lat',levels=levels, col='xofyear', col_wrap=3, extend='both')
            else:
                g=plot_data[mnth*6 : mnth*6+6,9,:,:].plot.contourf(x='lon', y='lat',levels=levels, col='xofyear', col_wrap=3, extend='both')
            plt.ylim(-30,60)
            plt.xlim(60,180)
            for i, ax in enumerate(g.axes.flat):
                land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False,ax=ax)
                #ax.quiver(data.lon[::3], data.lat[::1], data.u_av[mnth*6 +i,::1,::3], data.v_av[mnth*6+i,::1,::3], scale=qscale,headwidth=5)
                ax.set_title(pd_dic[mnth*6+i+1])
            plt.savefig('/scratch/rg419/plots/mom_budg_stat/'+ inp_fol + '/' + var + str(mnth + 1) + '.png')
            plt.close()
        
    vars = ['fv','dphidx', 'mom_trans', 'mom_stat','mom_sum']
    for var in vars:
        plot_mom_var(var,np.arange(-2.,2.1,0.2),500.)            


plot_pd_mom_fn('flat_10m')

