"""Make a seasonal mean of a field hemispherically symmetric"""
# 1/05/2018 Add function to flip field in latitude so that southern hemisphere and northern hemisphere switch places

import numpy as np
import xarray as xr

    
def make_sym(field, time_name='xofyear', asym=False):
    
    dims_in = list(field.dims)
    
    dims_extra = list(field.dims)    
    dims_extra.remove('lat')
    dims_extra.remove(time_name)

    
    field = field.transpose(time_name, 'lat', *dims_extra)

    # Convert to array
    field_values = field.values
    
    n = len(field.xofyear.values)//2
    field_out = np.zeros(field.values.shape)
    
    if asym:
        for i in range(0,n):
            field_out[i,:,...] = (field_values[i,:,...] - field_values[i+n,::-1,...])/2.
            field_out[i+n,:,...] = -1.*field_out[i,::-1,...]
    else:
        for i in range(0,n):
            field_out[i,:,...] = (field_values[i,:,...] + field_values[i+n,::-1,...])/2.
            field_out[i+n,:,...] = field_out[i,::-1,...]
    
        
    # Return to xarray
    field_out = xr.DataArray( field_out, dims = field.dims, coords = field.coords )
    field_out = field_out.transpose(*dims_in)
        
    return field_out


def flip_field(field, time_name='xofyear', asym=False):
    
    dims_in = list(field.dims)
    
    dims_extra = list(field.dims)    
    dims_extra.remove('lat')
    dims_extra.remove(time_name)

    
    field = field.transpose(time_name, 'lat', *dims_extra)

    # Convert to array
    field_values = field.values
    
    n = len(field.xofyear.values)//2
    field_out = np.zeros(field.values.shape)
    
    if asym:
        for i in range(0,n):
            field_out[i,:,...] = -1.*field_values[i+n,::-1,...]
            field_out[i+n,:,...] = -1.*field_values[i,::-1,...]
    else:
        for i in range(0,n):
            field_out[i,:,...] = field_values[i+n,::-1,...]
            field_out[i+n,:,...] = field_values[i,::-1,...]
    
        
    # Return to xarray
    field_out = xr.DataArray( field_out, dims = field.dims, coords = field.coords )
    field_out = field_out.transpose(*dims_in)
        
    return field_out
    

if __name__ == '__main__':
    """Examples/sanity check"""
    import matplotlib.pyplot as plt
    from hadley_cell import mass_streamfunction
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_1.000.nc')
    
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    
    precip_test = flip_field(data.precipitation)
    plt.figure(1)
    data.precipitation.mean('lon').plot.contourf(x='xofyear',y='lat')
    plt.figure(2)
    precip_test.mean('lon').plot.contourf(x='xofyear',y='lat')
    plt.show()
    
    precip_test = make_sym(data.precipitation)
    print('precip done')
    
    psi_test = make_sym(psi, asym=True)
    print('psi done')
    
    ucomp_test = make_sym(data.ucomp)
    print('u done')
    
    vcomp_test = make_sym(data.vcomp, asym=True)
    print('v done')
    
    plt.figure(1)
    precip_test.mean('lon').plot.contourf(x='xofyear',y='lat')
    
    plt.figure(2)
    psi_test.sel(pfull=500).plot.contourf(x='xofyear',y='lat')
    
    plt.figure(3)
    ucomp_test.mean('lon').sel(pfull=500).plot.contourf(x='xofyear',y='lat')
    
    plt.figure(4)
    vcomp_test.mean('lon').sel(pfull=900).plot.contourf(x='xofyear',y='lat')
    
    plt.show()