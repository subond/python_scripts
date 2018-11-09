"""Calculate gradients for scalars, vectors, products of vectors
   Assumes dimensions lat and lon in degrees, pfull in hPa
04/01/2018 - update so lon/lat/pfull dimension names can be specified as input - defaults remain the same"""

import numpy as np
import xarray as xr
from finite_difference import cfd#_old as cfd

def ddx(field, a = 6376.0e3, latname='lat', lonname='lon'):
    """Calculate d/dx of a given DataArray. DataArray must include a lat and lon dimensions"""
    
    try:
        field.coords[lonname]
    except:
        raise NameError('Coord ' + lonname + ' not found')
    try:
        field.coords[latname]
    except:
        raise NameError('Coord ' + latname + ' not found')
    
    coslat = np.cos(field[latname] * np.pi/180)
    field_dx = cfd( field.values, field[lonname].values*np.pi/180, field.get_axis_num(lonname) )   
    field_dx = xr.DataArray( field_dx, dims = field.dims, coords = field.coords )
    field_dx = field_dx/coslat/a
    
    return field_dx


def ddy(field, vector = True, uv = False, a = 6376.0e3, latname='lat'):
    """Calculate d/dy of a given DataArray. DataArray must include a lat dimension.
        kwargs: vector - specify if input field is vector or scalar
                prod   - if a vector, is the field uv?"""
    
    try:
        field.coords[latname]
    except:
        raise NameError('Coord ' + latname + ' not found')
        
    coslat = np.cos(field[latname] * np.pi/180)
    
    if vector and uv:
        cosfac = coslat**2
    elif vector:
        cosfac = coslat
    else:
        cosfac = 1.
    
    field_dy = cfd( (field*cosfac).values, field[latname].values*np.pi/180, field.get_axis_num(latname) )   
    field_dy = xr.DataArray( field_dy, dims = field.dims, coords = field.coords )
    field_dy = field_dy/cosfac/a
    
    return field_dy
    
    
def ddp(field, pname='pfull'):
    """Calculate d/dp of a given DataArray. DataArray must include a pfull dimension"""
    
    try:
        field.coords[pname]
    except:
        raise NameError('Coord ' + pname + ' not found')
    
    field_dp = cfd( field.values, field[pname].values*100., field.get_axis_num(pname) )   
    field_dp = xr.DataArray( field_dp, dims = field.dims, coords = field.coords )
    
    return field_dp


def ddt(field, timedir = 'xofyear', secperunit = 5.*86400., cyclic=True):
    """Calculate d/dt in unit/s of a given DataArray. DataArray must include a time dimension
       Define seconds per unit time using secperunit. Default calc is for pentads"""
    
    try:
        field.coords[timedir]
    except:
        raise NameError('Coord ' + timedir + ' not found')
    
    field_dt = cfd( field.values, field.coords[timedir].values*secperunit, field.get_axis_num(timedir), cyclic=cyclic)/secperunit/2.
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