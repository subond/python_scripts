#JP Streamfunction method
# 8/10/2018 Modify to calculate using divergent part of v if lons are specified, as this is the correct method. NB this may cause backwards compatibility issues where missing values are present at lower levels. For now, for aquaplanets in current use just ignore lowest 2 levels and add a switch to go back to original method. Will need to do production runs later. Will want to add vchi/uchi to datafiles as part of postprocessing in future, before interpolation

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def mass_streamfunction(data, a=6376.0e3, g=9.8, lons=[-1000], dp_in=0., intdown=True, use_v_locally=False):
    """Calculate the mass streamfunction for the atmosphere.

    Based on a vertical integral of the meridional wind.
    Ref: Physics of Climate, Peixoto & Oort, 1992.  p158.

    `a` is the radius of the planet (default Isca value 6376km).
    `g` is surface gravity (default Earth 9.8m/s^2).
    lons allows a local area to be used by specifying boundaries as e.g. [70,100]
    dp_in - if no phalf and if using regularly spaced pressure levels, use this increment for integral. Units hPa.
    intdown - choose integratation direction (i.e. from surface to TOA, or TOA to surface).

    Returns an xarray DataArray of mass streamfunction.
    """
    if lons[0]==-1000: #Use large negative value to use all data if no lons provided
        vbar = data.vcomp.mean('lon')
    elif use_v_locally:
        vbar = data.vcomp.sel(lon=lons).mean('lon')
    else:
        from windspharm.xarray import VectorWind
        # Create a VectorWind instance to handle the computation
        w = VectorWind(data.ucomp.sel(pfull=np.arange(50.,950.,50.)), data.vcomp.sel(pfull=np.arange(50.,950.,50.)))
        # Compute variables
        uchi, vbar, upsi, vpsi = w.helmholtz()
        vbar = vbar.sel(lon=lons).mean('lon').sortby('pfull', ascending=False)
        
    c = 2*np.pi*a*np.cos(vbar.lat*np.pi/180) / g
        
    # take a diff of half levels, and assign to pfull coordinates
    if dp_in==0.:
        dp=xr.DataArray(data.phalf.diff('phalf').values*100., coords=[('pfull', data.pfull)])
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
    #d = xarray.open_dataset('/scratch/jp492/gfdl_data/ref_earth/ref_earth_grey/run20/daily.nc',
    #    decode_times=False)
    #d['mass_sf'] = mass_streamfunction(d)
    #print((d.mass_sf))
    #d.close()
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/half_shallow.nc')
    lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= 170. and data.lon[i] < 190.]
    psi = mass_streamfunction(data, lons=lons, dp_in=50.)
    psi_old = mass_streamfunction(data, lons=lons, dp_in=50., use_v_locally=True)
    psi /= 1.e9
    psi_old /= 1.e9
    psi[:,40,:].plot.contour(x='lat', y='pfull', yincrease=False, levels=np.arange(0.,601,100.), colors='k', add_labels=False)
    psi[:,40,:].plot.contour( x='lat', y='pfull', yincrease=False, levels=np.arange(-600.,0.,100.), colors='k', linestyles='dashed', add_labels=False)
    plt.figure(2)
    psi_old[:,40,:].plot.contour(x='lat', y='pfull', yincrease=False, levels=np.arange(0.,601,100.), colors='k', add_labels=False)
    psi_old[:,40,:].plot.contour( x='lat', y='pfull', yincrease=False, levels=np.arange(-600.,0.,100.), colors='k', linestyles='dashed', add_labels=False)
    plt.show()