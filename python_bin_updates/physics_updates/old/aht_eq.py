"""
Calculate AHT_EQ, the atmospheric heat transport across the equator. Calculation based on Donohoe et al. 2013 (J. Clim.)

"""

import xarray as xr
from data_handling import cell_area
from physics import gradients as gr, model_constants as mc

def aht_eq(run):
    
    area = mc.a*mc.a*cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')   # Area of grid cells
    
    #Load in data, add area to dataset
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    data['area'] = (('lat','lon'), area)
    
    #Locate latitudes North and South of Equator
    lats_sh = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] <= 0]
    lats_nh = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= 0]
        
    # Calculate upward longwave flux at surface
    flux_lw_up = data.t_surf ** 4. * mc.stefan

    
    # Take global averages of:
    # Net SW at TOA +ve down
    toa_sw_anom = data.toa_sw - (data.toa_sw * data.area).sum(('lat','lon')) / data.area.sum(('lat','lon'))
    # Net SW at surface +ve down
    flux_sw_anom = data.flux_sw - (data.flux_sw * data.area).sum(('lat','lon')) / data.area.sum(('lat','lon'))
    # OLR +ve up
    olr_anom = data.olr - (data.olr * data.area).sum(('lat','lon')) / data.area.sum(('lat','lon'))
    # Net LW at surface +ve down
    flux_lw_anom = (data.flux_lw - flux_lw_up) - ( (data.flux_lw - flux_lw_up) * data.area).sum(('lat','lon')) / data.area.sum(('lat','lon'))
    # LH at surface +ve up
    flux_lhe_anom = data.flux_lhe - (data.flux_lhe * data.area).sum(('lat','lon')) / data.area.sum(('lat','lon'))
    # SENS at surface +ve up
    flux_t_anom = data.flux_t - (data.flux_t * data.area).sum(('lat','lon')) / data.area.sum(('lat','lon'))
    
    # Evaluate Atmospheric Heat Storage (AHS)
    ahs = gr.ddt( (data.temp * mc.cp_air + data.sphum * mc.L).sum('pfull')*5000./9.8 )
    ahs_anom = ahs - (ahs * data.area).sum(('lat','lon')) / data.area.sum(('lat','lon'))
    
    swabs = ((toa_sw_anom - flux_sw_anom) * data.area).sel(lat=lats_sh).sum(('lat','lon'))/1.e15
    olr = (olr_anom * data.area).sel(lat=lats_sh).sum(('lat','lon'))/1.e15
    shf = ((flux_t_anom + flux_lhe_anom - flux_lw_anom) * data.area).sel(lat=lats_sh).sum(('lat','lon'))/1.e15
    stor = (ahs_anom * data.area).sel(lat=lats_sh).sum(('lat','lon'))/1.e15
    
    ahteq = swabs - olr + shf - stor
        
    return ahteq, swabs, olr, shf, stor
    

