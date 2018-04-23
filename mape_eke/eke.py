# Load up climatology data and calculate EKE

import xarray as xr
from physics import model_constants as mc, gradients as gr
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as spint

def trop_height(temp, check=False):
    # Use temperature to calculate tropopause level
    # Take vertical gradient wrt height and find pressure level at which dT/dz goes and stays above -2 K/km
        
    pfull_big = np.arange(500.,0.,-10.)
    
    rho = temp.pfull * 100. / mc.rdgas / temp.mean('lon')
    dtdp = gr.ddp(temp.mean('lon'))
    dtdz = (-1. * dtdp * rho * mc.grav * 1000.)#.sel(pfull=levs)
    
    f = spint.interp1d(dtdz.pfull, dtdz,  axis=1, fill_value='extrapolate')
    dtdz_big = xr.DataArray(f(pfull_big), coords=[temp.xofyear, pfull_big, temp.lat], dims=['xofyear', 'pfull', 'lat'])
    
    trop = np.zeros([len(dtdz.xofyear), len(dtdz.lat)])
    for i in range(len(dtdz.xofyear)):
        for j in range(len(dtdz.lat)):
            trop[i,j] = min(pfull_big[np.where(dtdz_big[i,:,j] <= -2.)])
    
    trop = xr.DataArray(trop, coords=[temp.xofyear, temp.lat], dims=['xofyear', 'lat'])
    
    if check:
        trop.plot.contourf(levels=np.arange(90.,270.,10.))
        plt.figure(2)
        dtdz_big.sel(xofyear=68).plot.contourf(levels=np.arange(-12.,13.,1.))
        plt.show()
    
    return trop


def latmean(field, lats):
    coslat = np.cos(field.lat * np.pi/180.)
    field_mean = (field*coslat).sel(lat=lats).sum('lat') / coslat.sel(lat=lats).sum('lat')
    return field_mean
    

def eke(run):
    
    dp = 50. * 100.
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
        
    tke = (data.ucomp_sq + data.vcomp_sq).mean('lon')/2.
    zke = (data.ucomp.mean('lon')**2. + data.vcomp.mean('lon')**2.)/2.
    eke = tke - zke
    
    lowest_trop = trop_height(data.temp).max()
    #print 'Lowest tropopause level is ' + str(lowest_trop.values)
        
    strack_n = data.lat[data.lat >= 20.]
    strack_s = data.lat[data.lat <= -20.]
    levs = data.pfull[(data.pfull >= lowest_trop)]
    
    eke_n = latmean(eke, strack_n)
    eke_s = latmean(eke, strack_s)
        
    eke_n = eke_n.sel(pfull=levs).sum('pfull') * dp /mc.grav
    eke_s = eke_s.sel(pfull=levs).sum('pfull') * dp /mc.grav
    
    #print eke_n/1.e6
    #print eke_s/1.e6
    eke_out = (eke_n + eke_s)/2.
    return eke_n, eke_s, eke_out
    

#eke('full_qflux')
#eke('ap_20')
