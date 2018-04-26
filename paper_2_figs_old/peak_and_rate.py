"""
Find surface temperature centroid and rate of movement

"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_centroid
from physics import gradients as gr
from pylab import rcParams
from data_handling import cell_area
import scipy.interpolate as spint


    
def temp_max(data, lat_bound=45.):
    
    # Select latitudes over which to interpolate
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -lat_bound and data.lat[i] <= lat_bound]
    
    # Interpolate t_surf in latitude   
    f = spint.interp1d(lats, data.t_surf.sel(lat=lats).mean('lon').values, axis=-1, fill_value='extrapolate', kind='quadratic')
    
    lats_new = np.arange(-lat_bound, lat_bound+0.1, 0.1)
    t_new = f(lats_new)
    t_new = xr.DataArray(t_new, coords=[data.xofyear.values, lats_new], dims=['xofyear', 'lat'])
    
    locs = t_new.argmax('lat')

    t_max_lat = np.zeros(locs.values.shape)
    
    for i in range(len(locs)):
        t_max_lat[i] = t_new.lat[locs[i]].values
    
    t_max_lat = xr.DataArray(t_max_lat, coords=[data.xofyear.values], dims=['xofyear'])
    
    data['t_max_lat'] = t_max_lat
    
    return data
    
    


   
def precip_centroid_rate(run, ax_in, ylim_in=1.):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    precip_centroid(data)
    
    dpcentdt = gr.ddt(data.p_cent) * 86400.
    
    dpcentdt_ppos = dpcentdt.where(data.p_cent>=0.)
    dpcentdt_max = dpcentdt_ppos.where(dpcentdt_ppos==dpcentdt_ppos.max(),drop=True)
    pcent_dtmax = data.p_cent.sel(xofyear=dpcentdt_max.xofyear)
    print dpcentdt_max.values, pcent_dtmax.values   
    
    ax_twin = ax_in.twinx()
    
    data.p_cent.plot(ax=ax_twin, color='b')
    dpcentdt.plot(ax=ax_in, color='k')
        
    ax_in.set_xlabel('')
    ax_in.set_ylabel('')
    ax_twin.set_ylabel('')
    ax_twin.set_ylim([-30,30])
    ax_in.set_ylim([-1.*ylim_in, ylim_in])
    ax_in.set_title(run, fontsize=14)
    #ax_in.set_yticks([-8.,-4.,0.,4.,8.])
    ax_twin.spines['right'].set_color('blue')
    plt.tight_layout()
    
    return ax_twin



def pcent_grad_scatter(run, ax_in, ylim_in=1.):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    precip_centroid(data)
    dpcentdt = gr.ddt(data.p_cent) * 86400.
    
    ax_in.plot(data.p_cent, dpcentdt, 'xk')
        
    ax_in.set_xlabel('')
    ax_in.set_ylabel('')
    ax_in.set_xlim([-30,30])
    ax_in.set_ylim([-1.*ylim_in, ylim_in])
    ax_in.set_title(run, fontsize=14)
    
    

if __name__ == "__main__":
    
    plot_dir = '/scratch/rg419/plots/paper_2_figs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)

    rcParams['figure.figsize'] = 10, 12
    rcParams['font.size'] = 14
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')
    
    
    temp_max(data)
    data.t_max_lat.plot()
    plt.show()
    
    
    
    
    