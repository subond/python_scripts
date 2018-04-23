"""Calculate gradients for scalars, vectors, products of vectors
   Assumes dimensions lat and lon in degrees, pfull in hPa"""

import numpy as np
import xarray as xr
from finite_difference import cfd_old as cfd

def ddx(field, a = 6376.0e3):
    """Calculate d/dx of a given DataArray. DataArray must include lat and lon dimensions"""
    
    try:
        field.coords['lon']
    except:
        raise NameError('Coord lon not found')
    try:
        field.coords['lat']
    except:
        raise NameError('Coord lat not found')
    
    coslat = np.cos(field.lat * np.pi/180)
    field_dx = cfd( field.values, field.lon*np.pi/180, field.get_axis_num('lon') )   
    field_dx = xr.DataArray( field_dx, dims = field.dims, coords = field.coords )
    field_dx = field_dx/coslat/a
    
    return field_dx


def ddy(field, vector = True, uv = False, a = 6376.0e3):
    """Calculate d/dy of a given DataArray. DataArray must include lat dimension.
        kwargs: vector - specify if input field is vector or scalar
                prod   - if a vector, is the field uv?"""
    
    try:
        field.coords['lat']
    except:
        raise NameError('Coord lat not found')
        
    coslat = np.cos(field.lat * np.pi/180)
    
    if vector and uv:
        cosfac = coslat**2
    elif vector:
        cosfac = coslat
    else:
        cosfac = 1.
    
    field_dy = cfd( (field*cosfac).values, field.lat*np.pi/180, field.get_axis_num('lat') )   
    field_dy = xr.DataArray( field_dy, dims = field.dims, coords = field.coords )
    field_dy = field_dy/cosfac/a
    
    return field_dy
    
    
def ddp(field):
    """Calculate d/dp of a given DataArray. DataArray must include pfull dimension"""
    
    try:
        field.coords['pfull']
    except:
        raise NameError('Coord pfull not found')
    
    field_dp = cfd( field.values, field.pfull*100., field.get_axis_num('pfull') )   
    field_dp = xr.DataArray( field_dp, dims = field.dims, coords = field.coords )
    
    return field_dp


def ddt(field, timedir = 'xofyear', secperunit = 5.*86400.):
    """Calculate d/dt in unit/s of a given DataArray. DataArray must include a time dimension
       Define seconds per unit time using secperunit. Default calc is for pentads"""
    
    try:
        field.coords[timedir]
    except:
        raise NameError('Coord ' + timedir + ' not found')
    
    field_dt = cfd( field.values, field.coords[timedir].values*secperunit, field.get_axis_num(timedir) )   
    field_dt = xr.DataArray( field_dt, dims = field.dims, coords = field.coords )
    
    return field_dt



if __name__ == '__main__':
    """Examples/sanity check"""
    import matplotlib.pyplot as plt
    from data_handling import time_means
    
    data = time_means('full_qflux', [121,481], filename='plev_pentad', timeav='pentad')
    
    dudx = ddx(data.ucomp[40,17,:,:])
    dudy = ddy(data.ucomp[40,17,:,:])
    dudp = ddp(data.ucomp[40,:,:,:].mean('lon'))
    dudt = ddt(data.ucomp[:,17,32,0])
    
    data.ucomp[40,17,:,:].plot.contourf(levels=np.arange(-50.,51.,10.))
    plt.figure(2)
    dudx.plot.contourf(levels=np.arange(-0.000016,0.0000161,0.000002))
    plt.figure(3)
    dudy.plot.contourf(levels=np.arange(-0.000016,0.0000161,0.000002))
    
    plt.figure(4)
    data.ucomp[40,:,:,:].mean('lon').plot.contourf(levels=np.arange(-50.,51.,10.), yincrease=False)
    plt.figure(5)
    dudp.plot.contourf(levels=np.arange(-0.0045,0.0046,0.0005), yincrease=False)
    
    plt.figure(6)
    data.ucomp[:,17,32,0].plot()
    plt.figure(7)
    dudt.plot()
    
    plt.show()