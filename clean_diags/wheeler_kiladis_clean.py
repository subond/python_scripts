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
import scipy.signal

#Constants
rearth = 6376.0e3
omega = 7.2921150e-5
beta = 2.*omega/rearth  # in units [s^-1.m^-1]



def detrend_data(data, assym=False, region='tropics', period_fac=1., sanity_check=False):
    'Remove linear trend, then seasonal cycle from input data'
    if assym:
        data[:,0:len(data.lat)/2,:] = -1.*data[:,0:len(data.lat)/2,:]
        
    #Identify region to use
    lats = [data.lat[i] for i in range(len(data.lat)) if np.abs(data.lat[i]) <= 15]
    
    # Fit and remove any trend in annual mean
    data.coords['year'] = np.floor((data.time-1.)/(360.*period_fac))
    data_zav = data.mean('lon')
    data_annual = data_zav.groupby('year').mean(('time'))    
    data_local = data.sel(lat=lats)
        
    for i,l in enumerate(lats):
        p = np.polyfit(data_annual.year, data_annual.sel(lat=l).values, 1)        
        trend = (data_local.time-1.)/360. * p[0] + p[1]
        data_local[:,i,:] = data_local[:,i,:] - trend

    # Average the seasonal cycle and remove this
    data_local.coords['dayofyear'] = np.mod( data_local.time, 360.*period_fac)
    data_seasonal = data_local.groupby('dayofyear').mean(('time','lon'))
    data_seasonal = xr.DataArray( scipy.signal.savgol_filter(data_seasonal, 51, 3, axis=0), coords=data_seasonal.coords)

    data_detrended = data_local.groupby('dayofyear') - data_seasonal
    
    if sanity_check:
        plt.plot(data_local[0:360,0].mean('lon'))
        plt.plot(data_seasonal[:,0])
        plt.plot(data_detrended[0:360,0].mean('lon'))
        plt.show()
    
    return data_detrended


def get_JAS(data, numyears=30, dayrange=[180,270], sanity_check=False):
    'Convert time-lat-lon array to year-dayofyear-lat-lon, then select season using dayrange'
    
    data_reshaped = xr.DataArray( np.reshape(data.values,[numyears, 360, len(data.lat), len(data.lon)]),
                 [('year', np.arange(1.,numyears+1.)), ('dayofyear', data.dayofyear[0:360] ), ('lat', data.lat), ('lon', data.lon)]  )     
    
    if sanity_check:
        data[-1,:,:].plot.contourf()
        plt.figure(2)
        data_reshaped[29,359,:,:].plot.contourf()
        plt.show()
        
    return data_reshaped[:,dayrange[0]:dayrange[1]]


def taper(data, taperlen=10, sanity_check=False):
    'Taper the start and end of each season to prevent spectral leakage'
    
    if sanity_check:
        plt.plot(data[0,:,0,:].mean('lon'))
    
    data[:,:taperlen,:,:]  = data[:,:taperlen,:,:]  * (np.cos(np.linspace(-np.pi/2, 0, taperlen))**2)[np.newaxis, :, np.newaxis, np.newaxis]
    data[:,-taperlen:,:,:] = data[:,-taperlen:,:,:] * (np.cos(np.linspace(0, np.pi/2, taperlen))**2)[np.newaxis, :, np.newaxis, np.newaxis]
    
    if sanity_check:
        plt.plot(data[0,:,0,:].mean('lon'))
        plt.show()
    
    return data


def take_ffts(data, sanity_check=False, region='tropics', sym=None):
    
    lft = fft(data, axis=-1)     # FFT in space    
    tft = fft(lft, axis=1)               # FFT in time
    
    fts = np.fft.fftshift(tft, axes=(1,-1))
    # fourier transform in numpy is defined by exp(-2pi i (kx + wt)) 
    # but we want exp(kx - wt) so need to negate the x-domain
    fts = fts[:, :, :, ::-1]
    
    if sym is None:
        fts=fts
    elif sym is 'asym':
        fts = (fts - fts[:,:,::-1,:])/2.
    elif sym is 'sym':
        fts = (fts + fts[:,:,::-1,:])/2.
        
    power = np.abs(fts)**2  # Get power
    
    if region is 'tropics':
        spectra = power.mean(0).sum(1)  # Average power over all time windows, and sum over all latitudes
    elif region is 'ntropics':
        spectra = power.mean(0)[:,5:10,:].sum(1)  # NH only
    elif region is 'stropics':
        spectra = power.mean(0)[:,0:5,:].sum(1)  # SH only

    # Get wavenumber and frequency
    nw, nk = spectra.shape
    T = len(data.dayofyear)*1.  # Number of days used for spectra
    ks = np.arange(-nk/2, nk/2)  # non-dim wavenumber and frequency
    ws = np.arange(-nw/2, nw/2)
    ws = ws/T 
    
    spectra = xr.DataArray( spectra, [('ws', ws ), ('ks', ks )] )
    
    if sanity_check:
        levels = np.arange(14.,20.,0.2)
        np.log(spectra).plot.contourf(levels=levels)
        plt.xlim(-15.,15.)
        plt.ylim(0.,np.max(spectra.ws))
        plt.show()
    
    return spectra



def background(spectra, fsteps=10, ksteps=10, sanity_check=False):
    """Uses a 1-2-1 filter to generate 'red noise' background field for a spectra (as per WK1998)
        `fsteps` is the number of times to apply the filter in the frequency direction
        `ksteps` is the number of times to apply the filter in the wavenumber direction
    
    Returns a background field of same dimensions as `spectra`.
    """
    # create a 1D 1-2-1 averaging footprint
    bgf = spectra
    for i in range(fsteps):
        # repeated application of the 1-2-1 blur filter to the spectra
        footprint = np.array([[0,1,0], [0,2,0], [0,1,0]]) / 4.0
        bgf = scipy.signal.convolve2d(bgf, footprint, mode='same', boundary='wrap')
    for i in range(ksteps):
        # repeated application of the 1-2-1 blur filter to the spectra
        footprint = np.array([[0,0,0], [1,2,1], [0,0,0]]) / 4.0
        bgf = scipy.signal.convolve2d(bgf, footprint, mode='same', boundary='wrap')
    
    bgf = xr.DataArray( bgf, spectra.coords )
    
    if sanity_check:
        levels = np.arange(14.,20.,0.2)
        np.log(bgf).plot.contourf(levels=levels)
        plt.xlim(-15.,15.)
        plt.ylim(0.,np.max(spectra['cycles_per_day']))
        plt.show()
    
    return bgf


def remove_background(spectra, spectra_total=None):
    """A simple background removal to eliminate frequency noise.
       spectra: Symmetric or antisymmetric component of spectra
       spectra_total: Sum of symmetric and antisymmetric parts"""
    if spectra_total is None:
        spectra_total = spectra
    bg = background(spectra_total, fsteps=10, ksteps=10)
    return spectra / bg
    

def kelvin_dispersion(c, ks):
    k = ks/rearth
    w = c * k
    w = w * 86400. /2. /np.pi
    return w

def rossby_gravity_dispersion(c, ws, m):
    w = ws * 2. * np.pi / 86400.
    gkp = -(beta / (2*w)) + 0.5*np.sqrt((beta/w - 2*w/c)**2 - 8*m*beta/c)
    gkm = -(beta / (2*w)) - 0.5*np.sqrt((beta/w - 2*w/c)**2 - 8*m*beta/c)
    gkp = gkp * rearth
    gkm = gkm * rearth

    return gkp, gkm

def yanai_dispersion(c, ks):
    k = ks/rearth
    w = k*c/2 + 0.5*np.sqrt(k**2*c**2 + 4*beta*c)
    w = w * 86400. /2. /np.pi
    return w
    

def plot_wavelines(clist, mlist, ax=None, nk=False):
    inwk = np.linspace(-500., 500., 2000000)
    
    if ax is None:
        ax = plt.gca()
    
    for c in clist:
        if not nk:
            w_kelvin = kelvin_dispersion(c, inwk)
            ax.plot(inwk, w_kelvin, 'k')
        
        for m in mlist:
            if m != 0:
                gkp, gkm = rossby_gravity_dispersion(c, inwk, m)
                # Gravity waves: high frequency
                ax.plot(gkp[inwk > 0.2], inwk[inwk > 0.2], 'g')
                ax.plot(gkm[inwk > 0.2], inwk[inwk > 0.2], 'g')
                # Rossby waves: low frequency
                ax.plot(gkp[inwk < 0.2], inwk[inwk < 0.2], 'r')
                ax.plot(gkm[inwk < 0.2], inwk[inwk < 0.2], 'r')
            else:
                #Yanai Wave: when m = 0 only one solution is physically relevant
                w_yanai = yanai_dispersion(c, inwk)
                ax.plot(inwk, w_yanai, 'y')
                


def plot_wheeler_kiladis(run, clist, do_omega=False):
    rcParams['figure.figsize'] = 15, 6
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    if do_omega:
        #data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'_'+'omega500.0_daily.nc', decode_times=False)
        #v = detrend_data(data.omega)
        #figname = 'wk_omega_500.pdf'
        data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'_'+'ucomp150.0_daily.nc', decode_times=False)
        v = detrend_data(data.ucomp)
        figname = 'wk_u_150.pdf'
    else:
        data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'_'+'vcomp_daily.nc', decode_times=False)
        v = detrend_data(data.vcomp, assym=True)
        figname = 'wk_v_150.pdf'

    v_jas = get_JAS(v, dayrange=[180,270])

    v_taper = taper(v_jas)
        
    spectra = take_ffts(v_taper)
    
    spectra_asym = take_ffts(v_taper, sym='asym')
    spectra_sym = take_ffts(v_taper, sym='sym')
    
    spectra_nobg_a = remove_background(spectra_asym, spectra_total=(spectra_sym+spectra_asym)/2.)
    spectra_nobg_s = remove_background(spectra_sym, spectra_total=(spectra_sym+spectra_asym)/2.)
    
    #clist = np.sqrt(9.8*np.array([12.,25.,50.]))
    
    # Two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    levels = np.arange(1.,2.1,0.1)
    
    #First plot
    f1 = spectra_nobg_a.plot.contourf(
                 ax=ax1, levels = levels, cmap = 'inferno_r', add_colorbar=False, add_labels=False, extend='max')
    plot_wavelines(clist, [0, 2], ax=ax1, nk=True)
    ax1.set_xlim(-15.,15.)
    ax1.set_ylim(0.,np.max(spectra.ws))
    ax1.set_ylabel('Cycles per day')
    ax1.set_xlabel('Wavenumber')
    ax1.grid(True,linestyle=':')
    ax1.set_title('Antisymmetric modes')
    
    #Second plot
    f2 = spectra_nobg_s.plot.contourf(
                 ax=ax2, levels=levels, cmap = 'inferno_r', add_colorbar=False, add_labels=False, extend='max')
    plot_wavelines(clist, [1], ax=ax2)
    ax2.set_xlim(-15.,15.)
    ax1.set_ylim(0.,np.max(spectra.ws))
    ax2.set_xlabel('Wavenumber')
    ax2.grid(True,linestyle=':')
    ax2.set_title('Symmetric modes')    
    
    plt.subplots_adjust(right=0.95, top=0.9, bottom=0.15, hspace=0.1)
    plt.savefig(plot_dir+figname)
    plt.close()
    
plot_wheeler_kiladis('ap_2', [10.,50.]) 
plot_wheeler_kiladis('full_qflux', [10.,50.]) 
plot_wheeler_kiladis('ap_2', [10.,50.], do_omega=True) 
plot_wheeler_kiladis('full_qflux', [10.,50.], do_omega=True) 

