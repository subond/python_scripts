"""
Function to look at waveno spectra of v before and after onset at upper levels.

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from data_handling import month_dic
from pylab import rcParams
from numpy.fft import fft
        
def v_fft(run, months, before_after, filename='plev_daily', period_fac=1.):
    
    rcParams['figure.figsize'] = 10, 15
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True


    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Load in dataset
    name_temp = '/scratch/rg419/Data_moist/' + run + '/run%03d/'+filename+'.nc'
    names = [name_temp % m for m in range( months[0], months[1])  ]
    #read data into xarray 
    data = xr.open_mfdataset( names,
        decode_times=False,  # no calendar so tell netcdf lib
        # choose how data will be broken down into manageable chunks.
        chunks={'time': 30})
    data.coords['xofyear'] = np.mod( data.time, 360.*period_fac)
    
    vwnd = data.vcomp.load()

    # Take fourier transform and average over years
    N = 72*period_fac    # Number of samplepoints
    data.coords['wavenumber'] = np.arange(0,len(data.lon))    
    data['v_fft'] = (('time','lat','wavenumber'), 2.0/N * np.abs(fft(vwnd[:,17,:,:]))) #150 hPa
    data = data.groupby('xofyear').mean('time')    
    
    fft_before = data.v_fft[before_after[0]:before_after[0]+20,:,:16].mean('xofyear')
    fft_after = data.v_fft[before_after[1]:before_after[1]+20,:,:16].mean('xofyear')
    
    fft_5 = data.v_fft[:,33,:16]
    fft_45 = data.v_fft[:,48,:16]
    
    mn_dic = month_dic(1)
    tickspace = range(60,360,90)
    labels = [mn_dic[(k+30)/30 ] for k in tickspace]
    
    # Two subplots
    f, axarr = plt.subplots(2, sharex=True)
    plt.set_cmap('inferno_r')
    #First plot
    f1 = axarr[0].contourf(data.wavenumber[:16], data.lat, fft_before, levels = np.arange(3.,19.,3.), extend = 'max')
    axarr[0].set_ylabel('Latitude')
    axarr[0].set_yticks(np.arange(-60.,61.,30.))
    axarr[0].grid(True,linestyle=':')
    #Second plot
    f2 = axarr[1].contourf(data.wavenumber[:16], data.lat, fft_after, levels = np.arange(3.,19.,3.), extend = 'max')
    axarr[1].set_ylabel('Latitude')
    axarr[1].set_yticks(np.arange(-60.,61.,30.))
    axarr[1].set_xlabel('Wavenumber')
    axarr[1].grid(True,linestyle=':')
    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
    #Colorbar
    cb1=f.colorbar(f2, ax=axarr.ravel().tolist(), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
    cb1.set_label('Variance, $(m/s)^2$')
    
    plt.savefig(plot_dir+'fft_before_after.pdf', format='pdf')
    plt.close()
    
    
    
    rcParams['figure.figsize'] = 20, 10
    # Two subplots
    f, axarr = plt.subplots(2, sharex=True)
    #First plot
    f1 = fft_5.plot.contourf(x='xofyear', y='wavenumber', ax=axarr[0], levels = np.arange(3.,19.,3.), extend = 'max', add_colorbar=False, add_labels=False, cmap='inferno_r')
    axarr[0].set_ylabel('Wavenumber')
    axarr[0].grid(True,linestyle=':')
    #Second plot
    f2 = fft_45.plot.contourf(x='xofyear', y='wavenumber', ax=axarr[1], levels = np.arange(3.,19.,3.), extend = 'max', add_colorbar=False, add_labels=False, cmap='inferno_r')
    axarr[1].set_ylabel('Wavenumber')
    axarr[1].set_xticks(tickspace)
    axarr[1].set_xticklabels(labels,rotation=25)
    axarr[1].grid(True,linestyle=':')
    
    plt.subplots_adjust(right=0.95, top=0.95, bottom=0., hspace=0.1)
    #Colorbar
    cb1=f.colorbar(f2, ax=axarr.ravel().tolist(), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
    cb1.set_label('Variance, $(m/s)^2$')
    
    plt.savefig(plot_dir+'fft_5_45.pdf', format='pdf')
    plt.close()        
    
v_fft('ap_2', [121,481], [150,195])
v_fft('full_qflux', [121,481], [90,195])
v_fft('flat_qflux', [121,481], [90,220])


