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
    
rearth = 6376.0e3
omega = 7.2921150e-5
beta = 2.*omega/rearth  # in units [s^-1.m^-1]

def plot_wavelines(c, ks, ws, m=(1,3,5,7), kelvin=True, color=None, linestyle='--'):
    """Plot the analytic solutions to Kelvin, Rossby and gravity waves 
    on a wavenumber-frequency spectrum.
        `c` is wavespeed in [m.s^-1]
        `m` are the eigenvalues of the Rossby and gravity wave solutions to be plotted
    Plots lines in SI units"""

    k = ks / rearth                                 # in [wavelengths.m^-1]
    w = ws * 2. * np.pi / 86400.                                # in [wavelengths.s^-1]
    
    #khat = np.linspace(-500., 500., 2000000.)
    #what = np.linspace(-500., 500., 2000000.)
    
    #k = khat / rearth                                 # in [wavelengths.m^-1]
    #w = what * 2. * np.pi / 86400./5.                                # in [wavelengths.s^-1]
    
    if color is None:
        kcolor, gcolor, rcolor, ycolor = ('blue', 'black', 'green', 'red')
    else:
        try:
            kcolor, gcolor, rcolor, ycolor = color
        except:  # not a list of 4 colors, set all 4 to the same value
            kcolor = gcolor = rcolor = ycolor = color
    ls = linestyle
    
    # plot analytic solutions to equatorial waves
    # a. kelvin waves
    #   w = ck
    if kelvin:
        kline, = plt.plot(k, c*k, color=kcolor, linestyle=ls)
        kline.set_label('Kelvin Wave')

    # b. rossby waves
    # c. gravity waves
    #   w^2 - c^2 k^2 - \beta c^2 k / w = (2m+1) \beta c
    # The dispersion relation is quadratic in k and cubic in w
    # so solve for k and plot that way
    # plot several modes of m
    for mi in m:
        gkp = (-(beta / (2*w)) + 0.5*np.sqrt((beta/w - 2*w/c)**2 - 8*mi*beta/c)) 
        gkm = (-(beta / (2*w)) - 0.5*np.sqrt((beta/w - 2*w/c)**2 - 8*mi*beta/c)) 
        if mi != 0:
            # Gravity waves: high frequency
            gline, = plt.plot(gkp[w > 0.00001], w[w > 0.00001], color=gcolor, linestyle=ls)
            plt.plot(gkm[w > 0.00001], w[w > 0.00001], color=gcolor, linestyle=ls)
            # Rossby waves: low frequency
            rline, = plt.plot(gkp[w < 0.00001], w[w < 0.00001], color=rcolor, linestyle=ls)
            plt.plot(gkm[w < 0.00001], w[w < 0.00001], color=rcolor, linestyle=ls)
        else:
            # d. Yanai Wave: when m = 0 only one solution is physically relevant
            # w = kc/2 +- 1/2 sqrt(k^2 c^2 + 4 \beta c)
            yline, = plt.plot(k, k*c/2 + 0.5*np.sqrt(k**2*c**2 + 4*beta*c), color=ycolor, linestyle=ls)
            yline.set_label('Yanai Wave')

    gline.set_label('Gravity Waves')
    rline.set_label('Rossby Waves')
    
def axis_non_dim(c=10.0):
    non_dim_xticks = np.linspace(-10, 10, 11, dtype=np.int32)
    dim_xticks = non_dim_xticks / np.sqrt(c/beta)
    non_dim_yticks = np.linspace(0, 4, 5)
    dim_yticks = non_dim_yticks * np.sqrt(beta*c)
    plt.xticks(dim_xticks, non_dim_xticks)
    plt.yticks(dim_yticks, non_dim_yticks)
    plt.xlabel(r'Zonal Wavenumber $k \sqrt{c/\beta}$')
    plt.ylabel(r'Frequency $\omega / \sqrt{\beta c}$')

    

       
def kiladis_spectra(run, region='tropics', period_fac=1.):
    """Perform Wheeler-Kiladis Spectral Analysis on variable u.
        spinup: discard the first `spinup` days as initialisation
    Returns frequency-wavenumber spectra summed over tropical latitudes.
    """
    #Load in dataset (150 hPa meridional wind)
    v = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'_'+'vcomp_daily.nc', decode_times=False)
    
    #Identify region to sum over
    if region is 'tropics':
        lats = [v.lat[i] for i in range(len(v.lat)) if np.abs(v.lat[i]) <= 15]
    elif region is 'nmidlats':
        lats = [v.lat[i] for i in range(len(v.lat)) if v.lat[i] >= 30. and v.lat[i] <= 60.]
    elif region is 'smidlats':
        lats = [v.lat[i] for i in range(len(v.lat)) if v.lat[i] <= -30. and v.lat[i] >= -60.]
    elif region is 'ntropics':
        lats = [v.lat[i] for i in range(len(v.lat)) if v.lat[i] > 0 and v.lat[i] <= 15]
    elif region is 'stropics':
        lats = [v.lat[i] for i in range(len(v.lat)) if v.lat[i] < 0 and v.lat[i] >= -15]
    else:
        print 'Region not recognised, calculating over tropics'
        lats = [v.lat[i] for i in range(len(v.lat)) if np.abs(v.lat[i]) <= 15]
    
    v.vcomp[:,0:32,:] = -1.*v.vcomp[:,0:32,:]
    
    # Fit and remove any trend in annual mean
    v_zav = v.vcomp.mean('lon')
    v_zav.coords['year'] = np.floor((v.time-1.)/(360.*period_fac))
    v_annual = v_zav.groupby('year').mean(('time'))    
    v_local = v.vcomp.sel(lat=lats)
    
    for i,l in enumerate(lats):
        p = np.polyfit(v_annual.year, v_annual.sel(lat=l).values, 1)
        v_local[:,i] = v_local[:,i] - ( (v_local.time-1.)/360. * p[0] + p[1] )
    
    
    # Average the seasonal cycle and remove this
    v_local.coords['dayofyear'] = np.mod( v_local.time, 360.*period_fac)
    v_seasonal = v_local.groupby('dayofyear').mean(('time','lon'))
    v_detrended = v_local.groupby('dayofyear') - v_seasonal
    
    
    # Taper start and end of time series to prevent spectral leakage
    taper = 30
    v_detrended[:taper,:,:]  = v_detrended[:taper,:,:]  * (np.cos(np.linspace(-np.pi/2, 0, taper))**2)[:, np.newaxis, np.newaxis]
    v_detrended[-taper:,:,:] = v_detrended[-taper:,:,:] * (np.cos(np.linspace(0, np.pi/2, taper))**2)[:, np.newaxis, np.newaxis]
    
    #Take fourier transforms (copied from Jamie)
    lft = fft(v_detrended, axis=2)     # FFT in space
    tft = fft(lft, axis=0)               # FFT in time
    
    fts = np.fft.fftshift(tft, axes=(0,2))
    # fourier transform in numpy is defined by exp(-2pi i (kx + wt)) 
    # but we want exp(kx - wt) so need to negate the x-domain
    fts = fts[:, :, ::-1]
    # Sum over latitude
    fts = fts.sum(1)
    
    # Get wavenumber and frequency
    nw, nk = fts.shape
    T = (v.time[1] - v.time[0]).values*nw  #Get total length of time period: T = dt * no_vals
    ks = np.arange(-nk/2, nk/2)  # non-dim wavenumber and frequency
    ws = np.arange(-nw/2, nw/2)
    ws = ws/T
    #ksd = ks / rearth            # dimmed wavenumber and frequency (rad/sec)
    #wsd = ws * np.pi*2.
    
    fts = xr.DataArray( fts, [('cycles_per_day', ws ), ('wavenumber', ks )] )
    
    return fts
    
 
def power(spectra):
    return np.abs(spectra)**2
    

def background(spectra, fsteps=10, ksteps=10):
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
    
    return bgf

def remove_background(spectra, spectra_total=None):
    """A simple background removal to eliminate frequency noise.
       spectra: Symmetric or antisymmetric component of spectra
       spectra_total: Sum of symmetric and antisymmetric parts"""
    if spectra_total is None:
        spectra_total = spectra
    bg = background(spectra_total, fsteps=10, ksteps=10)
    return spectra - bg

spectra = kiladis_spectra('ap_2', region='tropics')
bg = background(spectra)
plt.contourf(bg)
plt.figure(2)
plt.contourf(spectra.values)
plt.show()

def plot_spectra_tropics(run):
    
    rcParams['figure.figsize'] = 15, 6
    rcParams['font.size'] = 25
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    spectra = kiladis_spectra(run, region='tropics')
    spectra_stropics = kiladis_spectra(run, region='stropics')
    spectra_ntropics = kiladis_spectra(run, region='ntropics')
    
    # Two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    levels = np.arange(11.,14.5,0.5)
    
    #First plot
    f1 = np.log(remove_background(spectra)).plot.contourf(
                 ax=ax1, levels = levels, cmap = 'viridis', add_colorbar=False, add_labels=False)
    ax1.set_xlim(-15.,15.)
    ax1.set_ylim(0.,np.max(spectra['cycles_per_day']))
    ax1.set_ylabel('Cycles per day')
    ax1.set_xlabel('Wavenumber')
    ax1.grid(True,linestyle=':')
    ax1.set_title('Symmetric modes')
    
    #Second plot
    f2 = np.log(remove_background(spectra_stropics - spectra_ntropics)).plot.contourf(
                 ax=ax2, levels=levels, cmap = 'viridis', add_colorbar=False, add_labels=False)
    ax2.set_xlim(-15.,15.)
    ax2.set_ylim(0.,np.max(spectra['cycles_per_day']))
    ax2.set_xlabel('Wavenumber')
    ax2.grid(True,linestyle=':')
    ax2.set_title('Antisymmetric modes')    
    
    plt.subplots_adjust(right=0.95, top=0.9, bottom=0.15, hspace=0.1)
    plt.show()
    
    plot_wavelines(10., spectra.wavenumber, spectra.cycles_per_day)
    axis_non_dim(c=10.)
    plt.xlim((-60/rearth, 60/rearth))    # show up to wavenumber +/- 100 
    plt.ylim((0, 2.0*np.pi/86400.))          # show positive frequencies, up to one wavelength per day
    plt.show()
    #plt.savefig(plot_dir + 'wheeler_kiladis.pdf', format='pdf')
    #plt.close()
        

def plot_spectra_midlats(run):
    #Sanity check: plot spectra for the NH midlatitudes    
    spectra_nmidlats = kiladis_spectra('full_qflux', region='nmidlats')
    
    np.log(remove_background(spectra_nmidlats)).plot.pcolormesh(levels=levels, cmap = 'viridis')
    plt.xlim(-15.,15.)
    plt.ylim(0.,np.max(spectra['Cycles per day']))
    plt.title('All modes (midlatitudes)')
    
    plt.show()
    

#plot_spectra_tropics('full_qflux')
plot_spectra_tropics('ap_2')


