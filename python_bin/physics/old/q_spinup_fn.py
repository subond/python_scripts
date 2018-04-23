#Calculate vertical integrals of area mean, annual mean specific humidity for whole atmosphere and levels above 100hPa (strat wv)
#Could set a threshold to determine spin-up end, e.g. require changes of less than 1% seems plausible. Need more data to confirm this is appropriate however

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import cell_area

def q_spinup(run, months, filename='atmos_pentad'):
    
    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    #read data into xarray 
    data = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    
    data.coords['year'] = data.time//360 + 1        
    data_yr = data.groupby('year').mean(('time'))
    
    area = cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')
    area_xr = xr.DataArray(area, [('lat', data.lat ), ('lon', data.lon)])
    dp = xr.DataArray( np.diff(data.phalf), [('pfull', data.pfull) ])
    
    p_strat = data.pfull[data.pfull <= 100.]
    
    #take area mean of q
    q_av = (data_yr.sphum*area_xr).sum(('lat','lon'))/area_xr.sum(('lat','lon'))

    #integrate over pressure levels above 100hPa and over whole atmosphere
   # q_strat = (q_avs[0:24,:]*dp[0:24]*100).sum(('pfull'))/9.8
    q_vint = (q_av*dp).sum('pfull')/9.8
    q_strat = (q_av*dp).sel(pfull=p_strat).sum('pfull')/9.8

    q_strat.plot()
    plt.xlabel('Year')
    plt.ylabel('Vertically integrated area mean specific humidity, kg/m^2')
    plotname = '/scratch/rg419/plots/spinup/qstrat_spinup_'+ run +'.png'
    plt.savefig(plotname)
    plt.close()
    
    q_vint.plot()
    plt.xlabel('Year')
    plt.ylabel('Vertically integrated area mean specific humidity, kg/m^2')
    plotname = '/scratch/rg419/plots/spinup/q_spinup_'+ run +'.png'
    plt.savefig(plotname)
    plt.close()

    return q_strat, q_vint



if __name__ == "__main__":

	q_strat, q_vint = q_spinup('ss_91.000', [61,250])
