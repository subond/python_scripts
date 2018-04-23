#load in a variable in pressure coordinates, plus theta
#Interpolate both to higher resolution
#Create bins
#Reorder columns

#For now, try just one snapshot

from physics import model_constants as mc
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def x_on_theta(run, field_list=['vcomp']):
    
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/plev_daily.nc'
    names = [name_temp % m for m in range( 121, 122)  ]
    
    data = xr.open_mfdataset( names, decode_times=False, chunks={'time': 30})
    
    #rC_fine = 98000:-1000:2000;
    #theta_grid = 290:5:350;
    
    theta = data.temp * (1000./data.pfull)**mc.kappa
    
    theta.mean(('lon','time')).plot.contourf()
    plt.show()
    
# Have not saved daily temperature so can not yet look at this. Bugger.

x_on_theta('ss_90.000')