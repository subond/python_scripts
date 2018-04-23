"""Add streamfunction, velocity potential, upsi, vpsi to file
   uchi and vchi can be calculated from u - upsi etc."""

import numpy as np
import os
from netCDF4 import Dataset
import xarray as xr
from windspharm.xarray import VectorWind


def add_vars_netcdf(filename, sf=True, vp=True, ups=True, vps=True):
    
    #first load up the dataset as an xarray and calculate the streamfunction
    data = xr.open_dataset(filename,decode_times=False)
    uwnd = data.ucomp
    vwnd = data.vcomp
    # Create a VectorWind instance to handle the computation
    w = VectorWind(uwnd, vwnd)

    # Compute variables
    streamfun, vel_pot = w.sfvp()
    uchi, vchi, upsi, vpsi = w.helmholtz()
    
    # Reverse latitude axes
    streamfun_out = streamfun.data[:,:,::-1,:]
    vel_pot_out = vel_pot.data[:,:,::-1,:]
    uchi_out = uchi.data[:,:,::-1,:]
    upsi_out = upsi.data[:,:,::-1,:]
    vchi_out = vchi.data[:,:,::-1,:]
    vpsi_out = vpsi.data[:,:,::-1,:]
    
    # Open file using netcdf
    dsin= Dataset(filename, 'a', format='NETCDF3_CLASSIC')
    
    #Add variables, updating if they already exist
    if sf:
        try:
            dsin['streamfun'][:] = streamfun_out
        except:
            sf_out = dsin.createVariable('streamfun', 'f4', ('time', 'pfull', 'lat', 'lon',))
            sf_out.setncatts({k: dsin.variables['ps'].getncattr(k) for k in dsin.variables['ps'].ncattrs()})
            sf_out.setncattr('long_name', 'streamfunction')
            sf_out.setncattr('units', 'm**2 s**-1')
            sf_out[:] = streamfun_out
    
    if vp:
        try:
            dsin['vel_pot'][:] = vel_pot_out
        except:
            vp_out = dsin.createVariable('vel_pot', 'f4', ('time', 'pfull', 'lat', 'lon',))
            vp_out.setncatts({k: dsin.variables['ps'].getncattr(k) for k in dsin.variables['ps'].ncattrs()})
            vp_out.setncattr('long_name', 'velocity potential')
            vp_out.setncattr('units', 'm**2 s**-1')
            vp_out[:] = vel_pot_out
    
    if ups:
        try:
            dsin['upsi'][:] = upsi_out
        except:
            sf_out = dsin.createVariable('upsi', 'f4', ('time', 'pfull', 'lat', 'lon',))
            sf_out.setncatts({k: dsin.variables['ps'].getncattr(k) for k in dsin.variables['ps'].ncattrs()})
            sf_out.setncattr('long_name', 'nondivergent zonal wind')
            sf_out.setncattr('units', 'm/sec')
            sf_out[:] = upsi_out
    
    if vps:
        try:
            dsin['vpsi'][:] = vpsi_out
        except:
            sf_out = dsin.createVariable('vpsi', 'f4', ('time', 'pfull', 'lat', 'lon',))
            sf_out.setncatts({k: dsin.variables['ps'].getncattr(k) for k in dsin.variables['ps'].ncattrs()})
            sf_out.setncattr('long_name', 'nondivergent meridional wind')
            sf_out.setncattr('units', 'm/sec')
            sf_out[:] = vpsi_out
        
    dsin.close()


if __name__ == "__main__":
    GFDL_DATA = os.environ['GFDL_DATA']
    filename_temp = GFDL_DATA + 'ap_1_rd/run%03d/atmos_pentad.nc'
    filenames = [filename_temp % m for m in range(1, 13)]
    
    for filename in filenames:
        print filename
        add_vars_netcdf(filename)

