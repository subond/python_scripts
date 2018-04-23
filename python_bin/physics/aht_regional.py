"""
Calculate the meridional atmospheric heat transport over a given region. 

"""

import xarray as xr
from data_handling import cell_area, rolling_mean
from physics import gradients as gr, model_constants as mc
import matplotlib.pyplot as plt
import numpy as np

def aht_regional(run, lonin=[-1.,361.], period_fac=1.):
    
    area = mc.a*mc.a*cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')   # Area of grid cells
    
    #Load in data, add area to dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['area'] = (('lat','lon'), area)
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    
    vh = data.vcomp_temp * mc.cp_air + data.sphum_v * mc.L + data.height*data.vcomp * mc.grav
    
    dvhdy = gr.ddy(vh)
    
    aht = ((dvhdy.sum('pfull')*5000./9.8) * data.area).sum(('lon')).cumsum('lat')
    
#    vh = ((gr.ddy(data.vcomp * h).sum('pfull')*5000./9.8)) #.sel(lon=lons).sum('lon')
    
   # aht_div = gr.ddy(vh)*data.area.
    
    #aht = aht_div.cumsum('lat')
    
    aht_rm = rolling_mean(aht, int(5*period_fac))
    
    aht_rm.plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.5e16,1.6e16,1.e15))
    plt.show()
    
    return aht_rm




if __name__ == '__main__':
    aht_regional('full_qflux_altevap2')

