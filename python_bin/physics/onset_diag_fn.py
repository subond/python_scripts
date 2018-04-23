#define monsoon onset as the pentad when rainfall of more than 8mm/day is observed
#read in precip data for a given run and locate onset pentad for each year at each gridpoint and produce pcolor plot

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import load_year_xr

def onset_diag_fn(rain):

    #mask values of precip less than 6mm/day
    tot_p_clim_masked = np.ma.masked_less(rain.values,8)

    #locate first unmasked value along pentad axis
    onset_index = np.ma.notmasked_edges(tot_p_clim_masked, axis=0)

    onset = np.zeros((64,128))
    onset[:] = np.nan
    #onset_end = np.zeros((64,128))
    #onset_end[:] = np.nan
    
    for i in range(0,len(onset_index[0][1])):
        onset[ onset_index[0][1][i], onset_index[0][2][i] ] = onset_index[0][0][i]+1
        #onset_end[ onset_index[1][1][i], onset_index[1][2][i] ] = onset_index[1][0][i]+1

    onset_pentad = xr.DataArray(onset, [('lat', rain.lat), ('lon', rain.lon)])
    #onset_end = xr.DataArray(onset_end, [('lat', rain.lat), ('lon', rain.lon)])

    return onset_pentad, onset_index#, onset_end

