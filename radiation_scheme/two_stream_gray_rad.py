"""Re-create two stream radiation code from GFDL-MiMA in offline python version"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.io import loadmat
from physics import model_constants as mc

# module constants
solar_constant  = 1360.0
del_sol         = 1.4
del_sw          = 0.0
ir_tau_eq       = 6.0
ir_tau_pole     = 1.5
atm_abs         = 0.0
sw_diff         = 0.0
linear_tau      = 0.1
wv_exponent     = 4.0
solar_exponent  = 4.0
albedo          = 0.3

rad_scheme = 'geen'

# constants for RG radiation version
ir_tau_co2_win  = 0.2150
ir_tau_wv_win1  = 147.11
ir_tau_wv_win2  = 1.0814e4
ir_tau_co2      = 1.5
ir_tau_wv1      = 41.0
ir_tau_wv2      = 104.
window          = 0.3732
carbon_conc     = 360.0

# parameters for Byrne and OGorman radiation scheme
bog_a = 0.8678
bog_b = 1997.9
bog_mu = 1.0


# Input profiles 
profiles = loadmat('q_and_t_sbdart.mat')

p_full = np.arange(2000.,100000.,4000.)
p_half = np.arange(0.,100001.,4000.)
yi = np.arange(-89.,90.,2.)
t = xr.DataArray(profiles['t_in'], [('lat', yi), ('pfull', p_full )])
q = xr.DataArray(profiles['q_in'], [('lat', yi), ('pfull', p_full )])

# Insolation
lat = t.lat * np.pi / 180.
p2          = (1. - 3.*np.sin(lat)**2)/4.
insolation  = 0.25 * solar_constant * (1.0 + del_sol * p2 + del_sw * np.sin(lat))


yy = t.lat * 90/60 * np.pi/180
sst = 27 * (1 - np.sin(yy)**2)
sst[abs(yi) >= 60] = 0

t_surf = sst + 273.15

b = mc.stefan*t**4
b_surf = mc.stefan*t_surf**4


n = len(t.pfull)



def sw_geen():
    # RG scheme: optical depth a function of wv and co2
    
    sw_down   = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    sw_up     = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    sw_tau_k  = xr.DataArray(np.zeros(t.lat.shape),       [('lat', yi)])
    sw_dtrans = xr.DataArray(np.zeros((len(t.lat), n)), [('lat', yi), ('pfull', p_full )])
    
        
    for k in range(0, n):
        sw_wv = sw_tau_k + 0.5194
        sw_wv = np.exp( 0.01887 / (sw_tau_k + 0.009522) + 1.603 / ( sw_wv*sw_wv ) )
        
        del_sol_tau = (( 0.0596 + 0.0029 * np.log(carbon_conc/360.)
                                + sw_wv[:] * q[:,k] ) 
                     * ( p_half[k+1] - p_half[k] ) / p_half[n])
        
        sw_dtrans[:,k] = np.exp( - del_sol_tau )
        sw_tau_k = sw_tau_k + del_sol_tau
            
    # compute downward shortwave flux
    sw_down[:,0] = insolation
    for k in range(0, n):
        sw_down[:,k+1]   = sw_down[:,k] * sw_dtrans[:,k]
    
    for k in range(0, n+1):
        sw_up[:,k]   = albedo * sw_down[:,n]
    
    return sw_down, sw_up
    

def sw_fb():
    #Frierson handling of SW radiation
    
    #SW optical thickness
    sw_tau_0    = (1.0 - sw_diff*np.sin(lat)**2)*atm_abs
    
    #compute optical depths for each model level
    sw_down   = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    sw_up     = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    sw_tau    = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    
    for k in range(0, n+1):
        sw_tau[:,k] = sw_tau_0 * (p_half[k]/mc.pstd_mks)**solar_exponent
    
    #compute downward shortwave flux
    for k in range(0, n+1):
        sw_down[:,k]   = insolation * np.exp(-sw_tau[:,k])
    
    for k in range(0, n+1):
        sw_up[:,k]   = albedo * sw_down[:,n]
        
    return sw_down, sw_up
    



def lw_down_geen():
    #split LW in 2 bands: water-vapour window and remaining = non-window
    #ref: Ruth Geen etal, GRL 2016 (supp. information].
    
    lw_dtrans     = xr.DataArray(np.zeros((len(t.lat), n)), [('lat', yi), ('pfull', p_full )])
    lw_dtrans_win = xr.DataArray(np.zeros((len(t.lat), n)), [('lat', yi), ('pfull', p_full )])
    lw_down       = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    lw_down_win   = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    
    for k in range(0, n):
        lw_del_tau    = (( ir_tau_co2 + 0.2023 * np.log(carbon_conc/360.)
                      + ir_tau_wv1 * np.log(ir_tau_wv2*q[:,k] + 1) )
               * ( p_half[k+1]-p_half[k] ) / p_half[n] )
        lw_dtrans[:,k] = np.exp( - lw_del_tau )
        
        lw_del_tau_win   = (( ir_tau_co2_win + 0.0954 * np.log(carbon_conc/360.)
                                     + ir_tau_wv_win1*q[:,k]
                                     + ir_tau_wv_win2*q[:,k]*q[:,k] )
                  * ( p_half[k+1]-p_half[k] ) / p_half[n] )
        lw_dtrans_win[:,k] = np.exp( - lw_del_tau_win )
        
    #compute downward longwave flux for window
    #Allocate a fraction of the longwave spectrum as window radiation
    b_win = window*b
    b_nw = (1.0 - window)*b
    #lw_down_win[:,0] = 0.0
    #lw_down[:,0] = 0.0
    
    for k in range(0, n):
        lw_down[:,k+1] = lw_down[:,k]*lw_dtrans[:,k] + b_nw[:,k]*(1. - lw_dtrans[:,k])
        
        lw_down_win[:,k+1] = (lw_down_win[:,k]*lw_dtrans_win[:,k]
                      + b_win[:,k]*(1.0 - lw_dtrans_win[:,k]))
    
    lw_down = lw_down + lw_down_win
        
    return lw_down, lw_dtrans, lw_dtrans_win


def lw_down_byrne():
    #dtau/ds = a*mu + b*q
    #ref: Byrne, M. P. & O'Gorman, P. A.
    #Land-ocean warming contrast over a wide range of climates:
    #Convective quasi-equilibrium theory and idealized simulations.
    #J. Climate 26, 4000-4106 (2013).
    
    lw_down   = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    lw_dtrans = xr.DataArray(np.zeros((len(t.lat), n)), [('lat', yi), ('pfull', p_full )])
    
    
    for k in range(0, n):
        lw_del_tau    = (bog_a*bog_mu + 0.17 * np.log(carbon_conc/360.)  + bog_b*q[:,k]) * (( p_half[k+1]-p_half[k] ) / p_half[n])
        lw_dtrans[:,k] = np.exp( - lw_del_tau )
        
    #compute downward longwave flux by integrating downward
    
    for k in range(0, n):
        lw_down[:,k+1] = lw_down[:,k]*lw_dtrans[:,k] + b[:,k]*(1. - lw_dtrans[:,k])
    
    return lw_down, lw_dtrans
    

def lw_down_frierson():
    #longwave optical thickness function of latitude and pressure
    
    lw_down   = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    lw_tau    = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    lw_dtrans = xr.DataArray(np.zeros((len(t.lat), n)), [('lat', yi), ('pfull', p_full )])
    
    lw_tau_0 = ir_tau_eq + (ir_tau_pole - ir_tau_eq)*np.sin(lat)**2
    
    #compute optical depths for each model level
    for k in range(0, n+1):
        lw_tau[:,k] = (lw_tau_0 * ( linear_tau * p_half[k]/mc.pstd_mks
                     + (1.0 - linear_tau) * (p_half[k]/mc.pstd_mks)**wv_exponent ))
    
    #longwave differential transmissivity
    for k in range(0, n):
        lw_dtrans[:,k] = np.exp( -(lw_tau[:,k+1] - lw_tau[:,k]) )
    
    #compute downward longwave flux by integrating downward
    lw_down[:,1]      = 0.
    for k in range(0, n):
        lw_down[:,k+1] = lw_down[:,k]*lw_dtrans[:,k] + b[:,k]*(1. - lw_dtrans[:,k])
    
    return lw_down, lw_dtrans
    

def lw_up_geen(lw_dtrans, lw_dtrans_win):
    #integrate upward, including window contribution
    
    lw_up     = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    lw_up_win = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    
    lw_up[:,n]     = b_surf*(1-window)
    lw_up_win[:,n] = b_surf*window
    
    b_win = window*b
    b_nw = (1.0 - window)*b
    
    for k in range(n-1,-1,-1):
        lw_up[:,k]   = lw_up[:,k+1]*lw_dtrans[:,k] + b_nw[:,k]*(1.0 - lw_dtrans[:,k])
        lw_up_win[:,k]   = lw_up_win[:,k+1]*lw_dtrans_win[:,k] + b_win[:,k]*(1.0 - lw_dtrans_win[:,k])
    
    lw_up = lw_up + lw_up_win
    return lw_up


def lw_up_fb(lw_dtrans):
    #compute upward longwave flux by integrating upward
    
    lw_up     = xr.DataArray(np.zeros((len(t.lat), n+1)), [('lat', yi), ('phalf', p_half )])
    
    lw_up[:,n]    = b_surf
    
    for k in range(n-1,-1,-1):
        lw_up[:,k]   = lw_up[:,k+1]*lw_dtrans[:,k] + b[:,k]*(1.0 - lw_dtrans[:,k])
    return lw_up    


if rad_scheme == 'frierson':
    sw_down, sw_up = sw_fb()
    lw_down, lw_dtrans = lw_down_frierson()
    lw_up = lw_up_fb(lw_dtrans)

elif rad_scheme == 'byrne':
    sw_down, sw_up = sw_fb()
    lw_down, lw_dtrans = lw_down_byrne()
    lw_up = lw_up_fb(lw_dtrans)    

elif rad_scheme == 'geen':
    sw_down, sw_up = sw_geen()
    lw_down, lw_dtrans, lw_dtrans_win = lw_down_geen()
    lw_up = lw_up_geen(lw_dtrans, lw_dtrans_win)    

else:
    print "Invalid scheme choice"

# net fluxes (positive up)
lw_flux  = lw_up - lw_down
sw_flux  = sw_up - sw_down
rad_flux = lw_flux + sw_flux

tdt_solar = xr.DataArray(np.zeros((len(t.lat), n)), [('lat', yi), ('pfull', p_full )])
tdt_rad   = xr.DataArray(np.zeros((len(t.lat), n)), [('lat', yi), ('pfull', p_full )])


for k in range(0, n):
    tdt_rad[:,k] = (( rad_flux[:,k+1] - rad_flux[:,k] )
        * mc.grav/( mc.cp_air*(p_half[k+1] - p_half[k]) ))
        
    tdt_solar[:,k] = (( sw_flux[:,k+1] - sw_flux[:,k] )
        * mc.grav/( mc.cp_air*(p_half[k+1] - p_half[k]) ))
    

cs = (tdt_solar*84600.).plot.contour(x='lat', y='pfull', yincrease = False, levels = np.arange(0,1.7,0.2), colors='k')
plt.clabel(cs, fmt='%1.1f')
plt.xlabel('Latitude')
plt.ylabel('Pressure, Pa')
plt.title('Shortwave heating')
plotname = '/scratch/rg419/plots/radiation_scheme/sw_htrt.png'
plt.savefig(plotname)
plt.close()
    
cs = ((tdt_rad - tdt_solar)*86400.).plot.contour(x='lat', y='pfull', yincrease = False, levels = np.arange(-10,10,0.5), colors='k')
plt.clabel(cs, fmt='%1.1f')
plt.xlabel('Latitude')
plt.ylabel('Pressure, Pa')
plt.title('Longwave heating')
plotname = '/scratch/rg419/plots/radiation_scheme/lw_htrt_150.png'
plt.savefig(plotname)
plt.close()
    