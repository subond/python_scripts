#Load in the radiative fluxes and look at how these varies with time - is there a signal from the strat wv?

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import cell_area

def ke_spinup(run, months, filename='atmos_pentad'):

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
    ke_av = ((data_yr.ucomp_sq + data_yr.vcomp_sq)*area_xr).sum(('lat','lon'))/area_xr.sum(('lat','lon'))

    #integrate over pressure levels above 100hPa and over whole atmosphere
    ke_vint = (ke_av*dp).sum('pfull')/9.8
    ke_strat = (ke_av*dp).sel(pfull=p_strat).sum('pfull')/9.8

    ke_strat.plot()
    plt.xlabel('Year')
    plt.ylabel('Vertically integrated area mean kinetic energy')
    plotname = '/scratch/rg419/plots/spinup/kestrat_spinup_'+ run +'.png'
    plt.savefig(plotname)
    plt.close()
    
    ke_vint.plot()
    plt.xlabel('Year')
    plt.ylabel('Vertically integrated area mean kinetic energy')
    plotname = '/scratch/rg419/plots/spinup/ke_spinup_'+ run +'.png'
    plt.savefig(plotname)
    plt.close()

    return ke_strat, ke_vint



if __name__ == "__main__":

	ke_strat, ke_vint = ke_spinup('full_qflux', [1,481])
	ke_strat, ke_vint = ke_spinup('ap_2', [1,481])
	ke_strat, ke_vint = ke_spinup('ap_20', [1,481])
    
    