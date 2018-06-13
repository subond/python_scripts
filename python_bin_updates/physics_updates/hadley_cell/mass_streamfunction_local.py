#Dubious local streamfunction method - basically mass_streamfunction but without the lon averaging. Remember this does not account for zonal motions

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

def mass_streamfunction_local(data, a=6376.0e3, g=9.8, lons=[-1000], dp_in=0., intdown=True):
    """Calculate the mass streamfunction for the atmosphere.

    Based on a vertical integral of the meridional wind.
    Ref: Physics of Climate, Peixoto & Oort, 1992.  p158.

    `a` is the radius of the planet (default Earth 6371km).
    `g` is surface gravity (default Earth 9.8m/s^2).

    Returns an xarray DataArray of mass streamfunction.
    """
    vbar = data.vcomp
    
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
    
    import sh
    from climatology import precip_centroid_ll

    data_dir = '/disca/restore/gv2scratch/rg419/Data_moist/climatologies/qflux_ss/'
    runs = ['qflux_0_100.nc', 'qflux_0_200.nc', 'qflux_0_300.nc',
            'qflux_5_100.nc', 'qflux_5_200.nc', 'qflux_5_300.nc',
            'qflux_10_100.nc', 'qflux_10_200.nc', 'qflux_10_300.nc',
            'qflux_15_100.nc', 'qflux_15_200.nc', 'qflux_15_300.nc',
            'qflux_20_100.nc', 'qflux_20_200.nc', 'qflux_20_300.nc',
            'qflux_25_100.nc', 'qflux_25_200.nc', 'qflux_25_300.nc']
    
    plot_dir = '/scratch/rg419/plots/qflux_ss/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    
    for run in runs:
        data = xr.open_dataset(data_dir + run)
        data['mass_sf'] = mass_streamfunction_local(data, dp_in=50.)
        data = precip_centroid_ll(data)
        (data.mass_sf/1.e9).sel(pfull=500.).plot.contourf(x='lon', y='lat', levels=np.arange(-700.,701.,100.))
        data.p_cent.plot.line('k')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        
        plt.savefig(plot_dir+'psi_500_' + run + '.pdf', format='pdf')
        plt.close()
        data.close()
        
        

    