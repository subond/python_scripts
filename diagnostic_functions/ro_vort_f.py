#load in vorticity output from model and produce seasonal plots of the zonal mean Rossby number vor/f

import numpy as np
import xarray as xr
from time_av_xr import load_syear_xr

#set run name
run_fol = 'aquaplanet_10m/np16'
inp_fol = 'aquaplanet_10m'
year = 5

rundata = load_syear_xr(run_fol, 5)
season_dic = {21.:'mam',
              22.:'jja',
              23.:'son',
              24.:'djf'}

rundata.coords['month'] = (rundata.time // 30) + 1
rundata.coords['season'] = ((rundata.time +30) // 90) 


#take time and zonal mean
vor_tzav = rundata.vor.groupby('season').mean(('lon', 'time'))
u_tzav = rundata.ucomp.groupby('season').mean(('lon', 'time'))

#set constants
omega = 7.2921150e-5
f = 2* omega * np.sin(rundata.lat *np.pi/180)

#evaluate rossby number
ro = -1*vor_tzav / f


#plot
ro.plot.contourf(x='lat', y='pfull',levels=np.arange(0,1.1,0.1), col='season', col_wrap=2, add_label = False, yincrease=False)

u_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-60,65,5), col='season', col_wrap=2, add_label = False, yincrease=False)

plt.show()
