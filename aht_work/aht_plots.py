"""
Calculate AHT_EQ, the atmospheric heat transport across the equator. Calculation based on Donohoe et al. 2013 (J. Clim.)
Plot this vs psi_eq at 500 hPa, precip lat, lat of 0 psi at 500 (i.e. reproduce Donohoe plots)
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
#import sh
from data_handling import cell_area
#from pylab import rcParams
from physics import aht_eq, model_constants as mc, mass_streamfunction
import scipy.interpolate as spint
from pylab import rcParams

month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 20
rcParams['text.usetex'] = True

def aht_vs_psi_eq(run, period_fac=1.):
    
    #Load in data, add area to dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')

    psi = mass_streamfunction(data, mc.a, dp_in=50.)
    psi /= 1.e9
    psi_eq = psi.sel(pfull=500.).isel(lat=[31,32]).mean('lat')
        
    ahteq = aht_eq(run)[0]
    
    ahteq.coords['month'] = np.mod( ahteq.xofyear -1, 72.*period_fac) //(6*period_fac) +1
    ahteq_month = ahteq.groupby('month').mean(('xofyear'))
    
    psi_eq.coords['month'] = np.mod( psi_eq.xofyear -1, 72.*period_fac) //(6*period_fac) +1    
    psi_eq_month = psi_eq.groupby('month').mean(('xofyear'))
    
    
    plt.plot(psi_eq, ahteq, 'xk', alpha=0.7)
    plt.plot(psi_eq_month, ahteq_month, 'xk', ms=7, mew=2)
    for i in range(0, len(month_list)):
        plt.text(psi_eq_month[i]+0.1, ahteq_month[i]+0.1, month_list[i], fontsize=14)
    plt.ylabel('Atmospheric heat transport at the Equator (PW)')
    plt.xlabel('500 hPa Mass Streamfunction at the Equator')
    plt.grid(True,linestyle=':')
    plt.xlim([-600,600])
    plt.ylim([-10,10])
    plt.savefig('/scratch/rg419/plots/aht_work/psi_aht_'+run+'.pdf', format='pdf')
    plt.close()


def precip_centroid(run, period_fac=1.):
    """
    Evaluate the precip centroid at each pentad. 
    """
    area = cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')   # Area of grid cells
    
    #Load in data, add area to dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['area'] = (('lat','lon'), area)
    
    # Get total precip
    try:
        data['precipitation'] = data.condensation_rain + data.convection_rain
    except:
        data['precipitation'] = data.precipitation
    
    # Select latitudes over which to evaluate precip centroid
    lat_bound = 20.
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -lat_bound and data.lat[i] <= lat_bound]
    
    # Integrate precip wrt longitude
    precip_area_lats = (data.precipitation.sel(lat=lats) * data.area.sel(lat=lats)).sum('lon').values
    
    # Interpolate precip in latitude    
    f = spint.interp1d(lats, precip_area_lats, axis=1, fill_value='extrapolate')
    lats_new = np.arange(-lat_bound, lat_bound+0.1, 0.1)
    p_new = f(lats_new)
    p_new = xr.DataArray(p_new, coords=[data.xofyear.values, lats_new], dims=['xofyear', 'lat'])
    
    # Calculate cumulative sum of precip with latitude
    p_area_int = p_new.cumsum('lat')
    
    # At each time find the precipitation centroid: the latitude at which half of the area integrated precip lies North/South
    p_cent = np.zeros((len(p_new.xofyear.values),))
    for i in range(1,len(p_new.xofyear.values)+1):
        p_cent[i-1] = p_new.lat[p_area_int.sel(xofyear=i) <= 0.5 * p_area_int.sel(xofyear=i).max('lat')].max('lat').values
        
    p_cent= xr.DataArray(p_cent, coords=[p_new.xofyear.values], dims=['xofyear'])
    
    # Calculate atmospheric heat transport at the equator
    ahteq = aht_eq(run)[0]
    
    # Calculate monthly averages
    ahteq.coords['month'] = np.mod( ahteq.xofyear -1, 72.*period_fac) //(6*period_fac) +1
    print ahteq.month
    ahteq_month = ahteq.groupby('month').mean(('xofyear'))
    
    p_cent.coords['month'] = np.mod( p_cent.xofyear -1, 72.*period_fac) //(6*period_fac) +1    
    p_cent_month = p_cent.groupby('month').mean(('xofyear'))
    
    
    # Plot
    plt.plot(ahteq, p_cent, 'xk', alpha=0.7)
    plt.plot(ahteq_month, p_cent_month, 'xk', ms=7, mew=2)
    for i in range(0, len(month_list)):
        plt.text(ahteq_month[i]+0.1, p_cent_month[i]+0.1, month_list[i], fontsize=14)
    plt.xlabel('Atmospheric heat transport at the Equator (PW)')
    plt.ylabel('Precipitation centroid latitude ($^{\circ}$)')
    plt.grid(True,linestyle=':')
    plt.ylim([-15,15])
    plt.xlim([-10,10])
    plt.savefig('/scratch/rg419/plots/aht_work/aht_pcent_'+run+'.pdf', format='pdf')
    plt.close()
    
   
aht_vs_psi_eq('full_qflux_altevap2')
aht_vs_psi_eq('sn_1.000')
aht_vs_psi_eq('rt_0.500')
aht_vs_psi_eq('rt_2.000')
aht_vs_psi_eq('sn_0.500', period_fac=0.5)
aht_vs_psi_eq('sn_2.000', period_fac=2.0)

#precip_centroid('full_qflux_altevap2')
#precip_centroid('sn_1.000')
#precip_centroid('rt_0.500')
#precip_centroid('rt_2.000')
#precip_centroid('sn_0.500', period_fac=0.5)
#precip_centroid('sn_2.000', period_fac=2.0)


