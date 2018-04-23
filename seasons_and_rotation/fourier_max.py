# fft the shit out of everything

from physics import mass_streamfunction
from data_handling import time_means
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from scipy.fftpack import fft


L = 2.500e6
cp =  287.04/(2./7.)
g = 9.8
stefan = 5.6734e-8

def fft_max(var, nmax=True, period_fac=1., plttitle=None):
    if nmax:
        var_max_loc = var.lat.values[np.argmax(var.mean('lon').values, axis=1)].tolist()
    else:
        var_max_loc = var.lat.values[np.argmin(var.mean('lon').values, axis=1)].tolist()
    
    # Number of samplepoints
    N = 72*period_fac
    yf = fft(var_max_loc)

    plt.plot( 2.0/N * np.abs(yf[:N/2]))
    plt.xlim(0,35)
    plt.xlabel('Wavenumber')
    plt.ylabel('Amplitude')
    plt.title(plttitle)
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+plttitle+'_fft.png')
    plt.close()
    
    plt.plot(var_max_loc)
    plt.xlabel('Pentad')
    plt.ylabel('Latitude')
    plt.title(plttitle)
    plt.savefig('/scratch/rg419/plots/seasons_and_rotation/'+plttitle+'_mean.png')
    plt.close()
        

data = time_means('sn_3.000',[121,481], filename='atmos_pentad', timeav='pentad', period_fac=3.)
fft_max(data.omega[:,27,:,:], nmax=False, period_fac=3., plttitle='sn_3.000_omega')
mse = (cp*data.temp + L*data.sphum + g*data.height)/1000.
fft_max(mse[:,38,:,:], period_fac=3., plttitle='sn_3.000_mse')
    
data = time_means('aquaplanet_10m',[121,481], filename='atmos_pentad', timeav='pentad')
fft_max(data.omega[:,27,:,:], nmax=False, plttitle='aquaplanet_10m_omega')
mse = (cp*data.temp + L*data.sphum + g*data.height)/1000.
fft_max(mse[:,38,:,:], plttitle='aquaplanet_10m_mse')

data = time_means('aquaplanet_2m',[121,481], filename='atmos_daily', timeav='pentad')
fft_max(data.omega[:,27,:,:], nmax=False, plttitle='aquaplanet_2m_omega')
mse = (cp*data.temp + L*data.sphum + g*data.height)/1000.
fft_max(mse[:,38,:,:], plttitle='aquaplanet_2m_mse')
