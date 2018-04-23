"""
Calculate AHT_EQ, the atmospheric heat transport across the equator. Calculation based on Donohoe et al. 2013 (J. Clim.)

"""

import xarray as xr
from data_handling import cell_area, rolling_mean
from physics import gradients as gr, model_constants as mc
import matplotlib.pyplot as plt
import numpy as np

def aht_zonal(run, period_fac=1.):
    
    area = mc.a*mc.a*cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')   # Area of grid cells
    
    #Load in data, add area to dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['area'] = (('lat','lon'), area)
        
    # Calculate upward longwave flux at surface
    flux_lw_up = data.t_surf ** 4. * mc.stefan
    
    # Evaluate Atmospheric Heat Storage (AHS)
    ahs = (gr.ddt( (data.temp * mc.cp_air + data.sphum * mc.L ).sum('pfull')*5000./9.8 ) * data.area).sum('lon')
    # SW absorbed by atmosphere
    swabs = ((data.toa_sw - data.flux_sw) * data.area).sum(('lon'))
    # Outgoing longwave
    olr = (data.olr * data.area).sum(('lon'))
    # Total upward surface heat flux
    shf = ((data.flux_t + data.flux_lhe - data.flux_lw + flux_lw_up) * data.area).sum(('lon'))
    
    aht_div = swabs - olr + shf - ahs
    aht = (aht_div - aht_div.mean('lat')).cumsum('lat')
    
    aht_rm = rolling_mean(aht, int(5*period_fac))
    
    aht_rm.plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.5e16,1.6e16,1.e15))
    plt.show()
    
    return aht_rm




if __name__ == '__main__':
    aht_zonal('sn_1.000')
    aht_zonal('sn_0.500')
    aht_zonal('sn_2.000')

