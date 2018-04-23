#JP Streamfunction method

import numpy as np
import xarray
import matplotlib.pyplot as plt

def mass_streamfunction(data, a=6376.0e3, g=9.8, lons=[-1000], dp_in=0., intdown=True):
    """Calculate the mass streamfunction for the atmosphere.

    Based on a vertical integral of the meridional wind.
    Ref: Physics of Climate, Peixoto & Oort, 1992.  p158.

    `a` is the radius of the planet (default Earth 6371km).
    `g` is surface gravity (default Earth 9.8m/s^2).

    Returns an xarray DataArray of mass streamfunction.
    """
    if lons[0]==-1000: #Use large negative value to use all data if no lons provided
        vbar = data.vcomp.mean('lon')
    else:
        vbar = data.vcomp.sel(lon=lons).mean('lon')
    c = 2*np.pi*a*np.cos(vbar.lat*np.pi/180) / g
        
    # take a diff of half levels, and assign to pfull coordinates
    if dp_in==0.:
        dp=xarray.DataArray(data.phalf.diff('phalf').values*100., coords=[('pfull', data.pfull)])
        return c*np.cumsum(vbar*dp, axis=vbar.dims.index('pfull'))
    else:
        dp=dp_in*100.
        if intdown:
            psi = c*dp* np.flip(vbar, axis=vbar.dims.index('pfull')) .cumsum('pfull')#[:,:,::-1]
        else:
            psi = -1.*c*dp* vbar .cumsum('pfull')#[:,:,::-1]
        #psi.plot.contourf(x='lat', y='pfull', yincrease=False)
        #plt.show()
        return psi #c*dp*vbar[:,::-1,:].cumsum('pfull')#[:,:,::-1]
        


if __name__ == '__main__':
    # example calculating EPV for a GFDL dataset
    d = xarray.open_dataset('/scratch/jp492/gfdl_data/ref_earth/ref_earth_grey/run20/daily.nc',
        decode_times=False)
    d['mass_sf'] = mass_streamfunction(d)
    print((d.mass_sf))
    d.close()