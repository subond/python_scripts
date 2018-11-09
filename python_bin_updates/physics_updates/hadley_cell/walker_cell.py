import numpy as np
import xarray
import matplotlib.pyplot as plt

def walker_cell(uchi, a=6376.0e3, latin=None, g=9.8, dp_in=0., intdown=False):
    """Calculate the Walker cell streamfunction for the atmosphere.

    Based on a vertical integral of the divergent component of the zonal wind.
    Ref: Nguyen et al. 2018 - Variability of the extent of the Hadley circulation in the southern hemisphere: a regional perspective
    
    uchi - divergent component of zonal wind
    `a` is the radius of the planet (default Isca value 6376km).
    `g` is surface gravity (default Earth 9.8m/s^2).
    lat - default is to look for the latitudes nearest the equator. lat allows input latitudes to be specified, over which an area weighted mean will be taken.
    dp_in - if no phalf and if using regularly spaced pressure levels, use this increment for integral. Units hPa.
    intdown - choose integratation direction along array

    Returns an xarray DataArray of mass streamfunction.
    """
    if latin == None:
        lats = uchi.lat.where(np.abs(uchi.lat)==np.abs(uchi.lat).min().values, drop=True)
    else:
        lats = [uchi.lat[i] for i in range(len(uchi.lat)) if uchi.lat[i] >= latin[0] and uchi.lat[i] <= latin[1]]
            
    coslat = np.cos(uchi.lat*np.pi/180.)
    uchi_areaweight = uchi * coslat
    ubar = uchi_areaweight.sel(lat=lats).sum('lat') / coslat.sel(lat=lats).sum('lat')
        
    c = 2*np.pi*a / g
        
    # take a diff of half levels, and assign to pfull coordinates
    if dp_in==0.:
        dp=xarray.DataArray(ubar.phalf.diff('phalf').values*100., coords=[('pfull', ubar.pfull)])
        return c*np.cumsum(ubar*dp, axis=ubar.dims.index('pfull'))
    else:
        dp=dp_in*100.
        if intdown:
            psi = c*dp* np.flip(ubar, axis=ubar.dims.index('pfull')) .cumsum('pfull')#[:,:,::-1]
        else:
            psi = -1.*c*dp* ubar.cumsum('pfull')#[:,:,::-1]
        #psi.plot.contourf(x='lat', y='pfull', yincrease=False)
        #plt.show()
        return psi #c*dp*vbar[:,::-1,:].cumsum('pfull')#[:,:,::-1]
        


if __name__ == '__main__':
    # example calculating Walker cell for a GFDL dataset
    from windspharm.xarray import VectorWind
    
    d = xarray.open_dataset('/disca/share/rg419/Data_moist/climatologies/3q_shallow.nc')
    
    w = VectorWind(d.ucomp.sel(pfull=np.arange(50.,950.,50.)), d.vcomp.sel(pfull=np.arange(50.,950.,50.)))
    # Compute variables
    uchi, vchi, upsi, vpsi = w.helmholtz()
    
    walker = walker_cell(uchi, dp_in=-50., intdown=True) 
    walker /= 1.e9
    walker.sel(xofyear=1).plot.contourf(x='lon', y='pfull', yincrease=False, levels=np.arange(-300.,301.,50.))
    plt.show()