# Fit optical depths using longwave fluxes

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import statsmodels.api as sm
import sh
from pylab import rcParams

# Load in one month's data
data = xr.open_dataset( '/scratch/rg419/Data_moist/rad_scheme/run120/atmos_daily.nc', decode_times=False)

# Choose plot folder    
plot_dir = '/scratch/rg419/plots/radiation_scheme/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)


# Define constants
STEFAN = 5.6734e-8
CP_AIR = 287.04 / (2./7.)


# Estimate fraction of blackbody radiation in the window using relation from thesis
window = (-0.0967* (data.temp/100.)**2. + 0.6516*(data.temp/100) - 0.7089).values
#window = ( data.uflx_win[:,40,...]/data.uflx[:,40,...] ).mean(('time', 'lat', 'lon')).values #use fixed fraction
#window = 0.37 #use constant used previously

# Calculate blackbody spectrum
bb = STEFAN * data.temp.values**4.
b_win = window*bb
b = (1.0 - window)*bb

# Evaluate estimate of non window optical depth from downward fluxes
dtau_nw = -np.log(( (data.dflx.values - data.dflx_win.values)[:,1:41, ...] - bb*(1.-window))/((data.dflx.values - data.dflx_win.values)[:,0:40, ...] - bb*(1.-window)))
print 'dtau_nw done'

# Evaluate estimate of window optical depth
dtau_win = -np.log((data.dflx_win.values[:,1:41, ...] - bb*window)/(data.dflx_win.values[:,0:40, ...] - bb*window))
print 'dtau_win done'

# Calculate layer thickness
data['dp'] =  (('pfull'), np.diff(data.phalf.values))	
dsigma = data.dp/data.ps*100.
dsigma = dsigma.transpose('time','pfull','lat','lon')

# Add optical depths to dataset
data['dtau_nw'] =  (('time','pfull','lat','lon'), dtau_nw/dsigma)	
data['dtau_win'] =  (('time','pfull','lat','lon'), dtau_win/dsigma)	
data['b_all'] =  (('time','pfull','lat','lon'), bb)	
data['b_nw'] =  (('time','pfull','lat','lon'), b)	
data['b_win'] =  (('time','pfull','lat','lon'), b_win)	


# Fit a ^0.5 relation to non window values
nw_exp = 0.5

q = data.sphum[:,10:40,...].values.flatten()
A = np.array([ q**nw_exp , np.ones(q.shape) ])
model = sm.OLS(data.dtau_nw[:,10:40,...].values.flatten(), A.T)
result=model.fit()
consts_nw = result.params
std_err_nw = result.bse

print consts_nw

# Fit a quadratic for window
A = np.array([ q**2, q , np.ones(q.shape) ])
model = sm.OLS(data.dtau_win[:,10:40,...].values.flatten(), A.T)
result=model.fit()
consts_win = result.params
std_err_win = result.bse

print consts_win


#Evaluate optical depths based on this fit
dtau_nw_s = consts_nw[0]*data.sphum**nw_exp + consts_nw[1]
dtau_win_s = consts_win[0]*data.sphum**2 + consts_win[1]*data.sphum + consts_win[2]

show_dtau_q = False
if show_dtau_q:
    # Plot versus humidity and look at relationship.
    rcParams['figure.figsize'] = 6, 7
    # Two subplots
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    #First plot
    ax1.plot(q, data.dtau_nw[:,10:40,...].values.flatten(), 'gx')
    ax1.plot(q, consts_nw[0]*q**0.5 + consts_nw[1], 'kx')
    ax1.set_ylabel('$d\tau/d\sigma$')
    ax1.grid(True,linestyle=':')

    ax2.plot(q, data.dtau_win[:,10:40,...].values.flatten(), 'gx')
    ax2.plot(q, consts_win[0]*q**2 + consts_win[1]*q + consts_win[2], 'kx')
    ax2.set_ylabel('$d\tau/d\sigma$')
    ax2.grid(True,linestyle=':')

    plt.subplots_adjust(left=0.17, right=0.83, top=0.95, bottom=0.1, hspace=0.2)
    plt.show()


# Plot lat-pressure optical depths to get idea of structure
rcParams['figure.figsize'] = 12, 7
# Four subplots
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

data.dtau_nw.mean(('lon','time')).plot.contourf(ax=ax1, levels=np.arange(0,12.1,1.), cmap='viridis', extend='both', add_labels=False)
ax1.invert_yaxis()
ax1.set_ylabel('Pressure')
ax1.set_title('$d\tau/d\sigma_{NW}$ RRTM')

data.dtau_win.mean(('lon','time')).plot.contourf(ax=ax2, levels=np.arange(0,12.1,1.), cmap='viridis', extend='both', add_labels=False)
ax2.invert_yaxis()
ax2.set_title('$d\tau/d\sigma_{WIN}$ RRTM')

dtau_nw_s.mean(('lon','time')).plot.contourf(ax=ax3, levels=np.arange(0,12.1,1.), cmap='viridis', extend='both', add_labels=False)
ax3.invert_yaxis()
ax3.set_ylabel('Pressure')
ax3.set_xlabel('Latitude')
ax3.set_title('$d\tau/d\sigma_{NW}$ GEEN')

dtau_win_s.mean(('lon','time')).plot.contourf(ax=ax4, levels=np.arange(0,12.1,1.), cmap='viridis', extend='both', add_labels=False)
ax4.invert_yaxis()
ax4.set_xlabel('Latitude')
ax4.set_title('$d\tau/d\sigma_{WIN}$ GEEN')

plt.savefig(plot_dir+'dtau_lat_p.pdf', format='pdf')
plt.close()
    

# Evaluate LW fluxes 
lw_down_win_s = np.zeros([30,41,64,128])
lw_down_s = np.zeros([30,41,64,128])
for k in range(0,40):
    lw_down_s[:,k+1,:,:] = lw_down_s[:,k,:,:]*np.exp( -1.*(dtau_nw_s*dsigma)[:,k,:,:]) + b[:,k,:,:]*(1. - np.exp( -1.*(dtau_nw_s*dsigma)[:,k,:,:]))
    lw_down_win_s[:,k+1,:,:] = lw_down_win_s[:,k,:,:]*np.exp( -1.*(dtau_win_s*dsigma)[:,k,:,:]) + b_win[:,k,:,:]*(1.0 - np.exp( -1.*(dtau_win_s*dsigma)[:,k,:,:]))

lw_up_win_s = np.zeros([30,41,64,128])
lw_up_s = np.zeros([30,41,64,128])
lw_up_s[:,40,:,:] = (data.uflx - data.uflx_win).values[:,40,:,:]
lw_up_win_s[:,40,:,:] = data.uflx_win.values[:,40,:,:]
for k in range(39,-1,-1):
    lw_up_s[:,k,:,:] = lw_up_s[:,k+1,:,:]*np.exp( -1.*(dtau_nw_s*dsigma)[:,k,:,:]) + b[:,k,:,:]*(1. - np.exp( -1.*(dtau_nw_s*dsigma)[:,k,:,:]))
    lw_up_win_s[:,k,:,:] = lw_up_win_s[:,k+1,:,:]*np.exp( -1.*(dtau_win_s*dsigma)[:,k,:,:]) + b_win[:,k,:,:]*(1.0 - np.exp( -1.*(dtau_win_s*dsigma)[:,k,:,:]))
    

# Evaluate heating rates
lw_flux  = lw_up_s + lw_up_win_s - lw_down_s - lw_down_win_s
lw_flux_rrtm  = data.uflx - data.dflx
lw_flux_win = lw_up_win_s -  lw_down_win_s
lw_flux_win_rrtm  = data.uflx_win - data.dflx_win

tdt_rad = np.zeros([30,40,64,128])
tdt_rad_rrtm = np.zeros([30,40,64,128])
tdt_rad_win = np.zeros([30,40,64,128])
tdt_rad_win_rrtm = np.zeros([30,40,64,128])
for k in range(0,40):
   tdt_rad[:,k,:,:]   = ( lw_flux[:,k+1,:,:] - lw_flux[:,k,:,:] ) * 9.8/( CP_AIR*data.dp[k].values*100. )
   tdt_rad_rrtm[:,k,:,:]   = ( lw_flux_rrtm[:,k+1,:,:] - lw_flux_rrtm[:,k,:,:] ) * 9.8/( CP_AIR*data.dp[k].values*100. )
   tdt_rad_win[:,k,:,:]   = ( lw_flux_win[:,k+1,:,:] - lw_flux_win[:,k,:,:] ) * 9.8/( CP_AIR*data.dp[k].values*100. )
   tdt_rad_win_rrtm[:,k,:,:]   = ( lw_flux_win_rrtm[:,k+1,:,:] - lw_flux_win_rrtm[:,k,:,:] ) * 9.8/( CP_AIR*data.dp[k].values*100. )


# Plot LW downward fluxes
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
im = ax1.contourf(data.lat, data.phalf, np.mean(np.mean(lw_down_s, axis=3), axis=0), levels = np.arange(0., 281., 40.), extend='both')
ax1.invert_yaxis()
ax1.set_ylabel('Pressure')
ax1.set_title('LWD$_{NW}$ GEEN')
plt.colorbar(im, ax=ax1)

im = ax2.contourf(data.lat, data.phalf, np.mean(np.mean(lw_down_win_s, axis=3), axis=0), levels = np.arange(0., 141., 20.), extend='both')
ax2.set_title('LWD$_{WIN}$ GEEN')
plt.colorbar(im, ax=ax2)


im = ax3.contourf(data.lat, data.phalf, np.mean(np.mean(data.dflx - data.dflx_win, axis=3), axis=0), levels = np.arange(0., 281., 40.), extend='both')
ax3.invert_yaxis()
ax3.set_xlabel('Latitude')
ax3.set_ylabel('Pressure')
ax3.set_title('LWD$_{NW}$ RRTM')
plt.colorbar(im, ax=ax3)

im = ax4.contourf(data.lat, data.phalf, np.mean(np.mean(data.dflx_win, axis=3), axis=0), levels = np.arange(0., 141., 20.), extend='both')
ax4.set_xlabel('Latitude')
ax4.set_title('LWD$_{WIN}$ RRTM')
plt.colorbar(im, ax=ax4)

plt.savefig(plot_dir+'lw_down.pdf', format='pdf')
plt.close()


# Plot LW upward fluxes
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
im = ax1.contourf(data.lat, data.phalf, np.mean(np.mean(lw_up_s, axis=3), axis=0), levels = np.arange(100., 275., 25.), extend='both')
ax1.invert_yaxis()
ax1.set_ylabel('Pressure')
ax1.set_title('LWU$_{NW}$ GEEN')
plt.colorbar(im, ax=ax1)

im = ax2.contourf(data.lat, data.phalf, np.mean(np.mean(lw_up_win_s, axis=3), axis=0), levels = np.arange(90., 210., 15.), extend='both')
ax2.set_title('LWU$_{WIN}$ GEEN')
plt.colorbar(im, ax=ax2)

im = ax3.contourf(data.lat, data.phalf, np.mean(np.mean(data.uflx - data.uflx_win, axis=3), axis=0), levels = np.arange(100., 275., 25.), extend='both')
ax3.invert_yaxis()
ax3.set_xlabel('Latitude')
ax3.set_ylabel('Pressure')
ax3.set_title('LWU$_{NW}$ RRTM')
plt.colorbar(im, ax=ax3)

im = ax4.contourf(data.lat, data.phalf, np.mean(np.mean(data.uflx_win, axis=3), axis=0), levels = np.arange(90., 210., 15.), extend='both')
ax4.set_xlabel('Latitude')
ax4.set_title('LWU$_{WIN}$ RRTM')
plt.colorbar(im, ax=ax4)

plt.savefig(plot_dir+'lw_up.pdf', format='pdf')
plt.close()


# Plot heat rates
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
im = ax1.contourf(data.lat, data.pfull, np.mean(np.mean(tdt_rad, axis=3), axis=0), levels = np.arange(-0.00012, 0.000021, 0.00001), extend='both')
ax1.invert_yaxis()
ax1.set_ylabel('Pressure')
ax1.set_title('Total heatrate GEEN')
plt.colorbar(im, ax=ax1)

im = ax2.contourf(data.lat, data.pfull, np.mean(np.mean(tdt_rad_rrtm, axis=3), axis=0), levels = np.arange(-0.00012, 0.000021, 0.00001), extend='both')
ax2.set_title('Total heatrate RRTM')
plt.colorbar(im, ax=ax2)

im = ax3.contourf(data.lat, data.pfull, np.mean(np.mean(tdt_rad_win, axis=3), axis=0), levels = np.arange(-0.000018, 0.0000041, 0.000002), extend='both')
ax3.invert_yaxis()
ax3.set_xlabel('Latitude')
ax3.set_ylabel('Pressure')
ax3.set_title('Window heatrate GEEN')
plt.colorbar(im, ax=ax3)

im = ax4.contourf(data.lat, data.pfull, np.mean(np.mean(tdt_rad_win_rrtm, axis=3), axis=0), levels = np.arange(-0.000018, 0.0000041, 0.000002), extend='both')
ax4.set_xlabel('Latitude')
ax4.set_title('WIindow heatrate RRTM')
plt.colorbar(im, ax=ax4)

plt.savefig(plot_dir+'htrt.pdf', format='pdf')
plt.close()