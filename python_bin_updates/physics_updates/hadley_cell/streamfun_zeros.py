"""
Find the latitudes at which the Hadley cell streamfunction changes sign at a given pressure 22/02/2018
Similar to psi_edge_loc, but for simpler cases with no missing values 
"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from hadley_cell import mass_streamfunction

# Define latitudes to interpolate psi onto
lats = np.arange(-60, 60.1, 0.5)


def get_streamfun_zeros(data, lev=150., sanity_check=False):
    # Calculate mass streamfunction
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    psi = psi.sel(pfull=lev)
    
    n = len(psi.xofyear.values) 
     
    psi_mask = np.ones(psi.shape)
    psi_mask[psi <= 0] = 0.
    
    psi_mask = xr.DataArray(psi_mask, coords=psi.coords, dims=['lat', 'xofyear'])
    
    borders = psi_mask.diff('lat')
    
    # Identify Hadley cell edges and ITCZ by:
    # ITCZ: positive border closest to the equator
    # cell edges: negative border poleward of ITCZ
    
    lat_itcz = np.zeros(n,)
    lat_nh = np.zeros(n,)
    lat_sh = np.zeros(n,)
    for i in range(n):
        lats_pos = data.latb[1:-1].values[np.where(borders[:,i] == 1)[0]]
        lats_pos = lats_pos[(lats_pos>-30.) & (lats_pos<30.)]
        
        if len(lats_pos) > 1:
            if np.all(lats_pos>0.):
                lats_pos = np.max(lats_pos)
            elif np.all(lats_pos<0.):
                lats_pos = np.min(lats_pos)
            else:
                lats_pos = lats_pos[ np.abs(lats_pos)==np.max(np.abs(lats_pos))]
        lat_itcz[i] = lats_pos
        
        lats_neg = data.latb[1:-1].values[np.where(borders[:,i] == -1)[0]]
        lat_nh[i] = np.min(lats_neg[(lats_neg - lat_itcz[i] > 0) & (lats_neg > 0)])
        lat_sh[i] = np.max(lats_neg[(lats_neg - lat_itcz[i] < 0) & (lats_neg < 0)])
        
    if sanity_check == True:
        psi.plot.contourf()
        plt.plot(data.xofyear, lat_itcz,'k')
        plt.plot(data.xofyear, lat_nh,'k')
        plt.plot(data.xofyear, lat_sh,'k')
        plt.show()
    
    return lat_itcz, lat_nh, lat_sh
    
    
    
    


#data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_0.500.nc')
#get_streamfun_zeros(data)

#data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')
#get_streamfun_zeros(data)

#data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/rt_2.000.nc')
#get_streamfun_zeros(data)