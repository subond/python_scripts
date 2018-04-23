# Load up climatology data and calculate MAPE

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
    

def mape(run):
    # NB results are larger than in O'Gorman and Schneider papers, but are similar to those in Li et al 2007 - Lorenz Energy Cycle of the Global Atmosphere (GRL), and Lambert 1984 - A global available potential energy-kinetic energy budget (Atmos-Ocean), Saltzman 1970 - Large-Scale Atmospheric Energetics in the Wave Number Domain (rev Geophys), Peixoto and Oort 1974 - The annual distribution of atmospheric energy on a planetary scale.
    # Get smaller values if confine average over a smaller storm track region and look at annual mean - might be what's going on in OG + S, but overall I think my code is fine. May get smaller values in my perpetual equinox runs for example
    
    p00 = 1000.
    dp = 50.
    prefac = mc.cp_air * p00 * 100. / 2. /mc.grav
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
        
    theta = data.temp * (p00 / data.pfull) ** mc.kappa
    
    lowest_trop = trop_height(data.temp).max()
    #print 'Lowest tropopause level is ' + str(lowest_trop.values)
    
    theta_zav = theta.mean('lon')
    thetasq_zav = theta_zav ** 2.
    dthdp = gr.ddp(theta_zav)
    
    strack_n = theta.lat[theta.lat >= 20.]
    strack_s = theta.lat[theta.lat <= -20.]
    levs = data.pfull[(data.pfull >= lowest_trop) & (data.pfull <= 900.)]
    
    theta_n = latmean(theta_zav, strack_n)
    theta_s = latmean(theta_zav, strack_s)
    thetasq_n = latmean(thetasq_zav, strack_n)
    thetasq_s = latmean(thetasq_zav, strack_s)
    dthdp_n = latmean(dthdp, strack_n)
    dthdp_s = latmean(dthdp, strack_s)
    
    gamma_n = -1. * mc.kappa / (data.pfull * 100.) / dthdp_n
    gamma_s = -1. * mc.kappa / (data.pfull * 100.) / dthdp_s
    
    theta_var_n = thetasq_n - theta_n**2.
    theta_var_s = thetasq_s - theta_s**2.
    
    mape_integrand_n = (data.pfull / p00) ** mc.kappa * gamma_n * theta_var_n
    mape_integrand_s = (data.pfull / p00) ** mc.kappa * gamma_s * theta_var_s
    
    #print mape_integrand_s[:,0]
    
    mape_n = prefac * mape_integrand_n.sel(pfull=levs).sum('pfull') * dp /p00
    mape_s = prefac * mape_integrand_s.sel(pfull=levs).sum('pfull') * dp /p00
    
    #print mape_n/1.e6
    #print mape_s/1.e6
    mape = (mape_n + mape_s)/2.
    return mape_n, mape_s, mape
    

#mape('full_qflux')
#mape('ap_20')
