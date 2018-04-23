"""
Load up climatology for an experiment including land and plot the mean SST over land and over ocean as a function of time.

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from pylab import rcParams
from stephen import sagp

def area_average(data, variable_name, model_params, latin=[-90.,90.], lonin=[-1.,361.], land_ocean_all='all', level=None, axis_in='xofyear'):
    
    print 'performing area average on ',variable_name, 'of type ', land_ocean_all
    
    if lonin[1]>lonin[0]:
        lons = data.lon[ np.array(data.lon >= lonin[0]) & np.array(data.lon < lonin[1]) ]
    else:
        lons = data.lon[ np.array(data.lon >= lonin[0]) | np.array(data.lon < lonin[1]) ]
        
    lats = data.lat[ np.array(data.lat >= latin[0]) & np.array(data.lat < latin[1]) ]
        
    data['region_mask'] = (('lat','lon'), np.zeros_like(data.land))
    
    data.region_mask.loc[dict(lon = lons, lat = lats)] = 1.
    
    if(level!=None):
        data_to_average=data[variable_name].sel(pfull=level, method='nearest')
    else:
        data_to_average=data[variable_name]
        
    try:
        grid_area=data['grid_cell_area']
    except KeyError:
        sagp.get_grid_sizes(data,model_params)
        grid_area=data['grid_cell_area']
        
    if(land_ocean_all == 'land'):
        scaled_grid_area=grid_area*(data['land'])*data.region_mask

    elif(land_ocean_all == 'ocean'):
        scaled_grid_area=grid_area*(1.-data['land'])*data.region_mask
        
    elif(land_ocean_all == 'all'):
        scaled_grid_area=grid_area*data.region_mask
    
    else:
        print 'invalid area-average option: ',land_ocean_all
        return
    
    #scaled_grid_area.plot.contourf()
    #plt.show()
    
    multiplied=scaled_grid_area*data_to_average
    average=multiplied.sum(('lat','lon'))/scaled_grid_area.sum(('lat','lon'))
        
    new_var_name=variable_name+'_area_av_'+land_ocean_all
    
    data[new_var_name]=((axis_in), average)
    
	


def sst_land_ocean(run, lonin=[-1.,361.], latin=[-90.,90.]):
    
    #rcParams['figure.figsize'] = 10, 7
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')
    data['land'] = (('lat','lon'),land.land_mask.values)
    
    model_params = sagp.model_params_set('/scratch/rg419/GFDL_model/GFDLmoistModel/', delta_t=720., ml_depth=20.)
        
    area_average(data, 't_surf', model_params, land_ocean_all='land', axis_in='xofyear', lonin=lonin, latin=latin)
    area_average(data, 't_surf', model_params, land_ocean_all='ocean', axis_in='xofyear', lonin=lonin, latin=latin)
    
    data.t_surf_area_av_land.plot()
    data.t_surf_area_av_ocean.plot()
    plt.show()
    

sst_land_ocean('full_qflux', latin=[0.,30.], lonin=[60.,150.])
sst_land_ocean('flat_qflux', latin=[0.,30.], lonin=[60.,150.])
