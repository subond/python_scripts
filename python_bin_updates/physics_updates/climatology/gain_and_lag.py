# 13/12/2017 Calculate lag and gain of SST relative to insolation

from data_handling_updates import month_dic, gradients as gr, model_constants as mc
import numpy as np
from numpy.fft import fft
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh


rcParams['figure.figsize'] = 6, 10
rcParams['font.size'] = 16
plot_dir = '/scratch/rg419/plots/surface_fluxes/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


def fft_test(fn):
    # sanity check behaviours of fft functions
    
    test = fft(fn)
    test_amp = np.abs(test) * 2./len(times)
    test_phase = np.arctan(np.imag(test)/np.real(test)) * 180./np.pi

    plt.figure(1)
    plt.plot(fn)
    plt.figure(2)
    plt.plot(test_amp[0:20])
    plt.figure(3)
    plt.plot(test_phase[0:20])
    plt.show()
    

def get_fft(data, sanity_check=False):
    #Take a dataset, remove the time mean (dim 0) and take an fft. Return amplitude and phase
    
    data_anom = data - data.mean('xofyear')
    
    data_fft = fft(data_anom, axis=0)
    
    data_amp = np.abs(data_fft) * 2. / data.shape[0]
    
    data_phase = np.arctan(np.imag(data_fft)/np.real(data_fft)) * 180./np.pi
    
    if sanity_check:
    
        plt.figure(1)
        data_anom.mean('lon').plot.contourf(x='xofyear', y='lat')    
    
        plt.figure(2)
        plt.contourf(np.mean(data_amp[0:20,:,:],2).T)
        plt.colorbar()
    
        plt.figure(3)
        plt.contourf(np.mean(data_phase[0:20,:,:],2).T)
        plt.colorbar()
        
        plt.show()
    
    return data_amp, data_phase
    

def gain_and_lag(run):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    sst_amp, sst_phase = get_fft(data.t_surf)
    sw_amp, sw_phase = get_fft(data.flux_sw)
    
    gain = sst_amp[1,:,:]/sw_amp[1,:,:]
    lag = -1.*(sst_phase[1,:,:] - sw_phase[1,:,:])
    
    gain = xr.DataArray(gain, coords=[data.lat, data.lon], dims=['lat', 'lon'])
    lag = xr.DataArray(lag, coords=[data.lat, data.lon], dims=['lat', 'lon'])
    
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    
    gain.plot.contourf(ax=ax1, x='lon', y='lat', levels=np.arange(0.,0.21,0.01), extend='neither', add_labels=False, cmap='RdBu_r')
    
    lag.plot.contourf(ax=ax2, x='lon', y='lat', levels=np.arange(0.,91.,5.), extend='neither', add_labels=False, cmap='RdBu_r')
    
    plt.savefig(plot_dir + 'gain_and_lag_' + run + '.pdf', format='pdf')
    plt.close()    
    

if __name__ == "__main__":


    gain_and_lag('half_shallow')
    gain_and_lag('half_shallow_5')
    gain_and_lag('half_shallow_10')

    
    
    
    #times = np.arange(0,360.)*np.pi/180.
    #fft_test(np.sin(times))
