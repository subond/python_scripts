#Load in u, v, omega, geopotential, evaluate the terms in the zonal momentum budget
#NB may also need to account for damping terms

#Start by doing this for an annual average, can extend after

import sys
sys.path.append('/scratch/rg419/workdir_moist/python_bin')
from time_av_xr import load_year_xr
import numpy as np
import matplotlib.pyplot as plt
from StreamFunction import cal_stream_fn
import xarray as xr

#load in run data
run_fol = 'aquaplanet_10m/np16'
year = 7
rundata = load_year_xr(run_fol, [year])

data = {'omega':np.ma.asarray((rundata.omega.mean(('lon'))).values),
        'vcomp':np.ma.asarray((rundata.vcomp.mean(('lon'))).values),
        'lat':rundata.lat.values,
        'time':rundata.time.values,
        'pfull':rundata.pfull.values}

psi = cal_stream_fn(data,1)
psi.shape

plt.figure(1)
plt.contourf(rundata.lat,rundata.pfull,np.mean(psi,0))
plt.colorbar()

plt.figure(2)
plt.plot(rundata.lat,np.mean(psi[:,0,:],0))

plt.figure(3)
plt.plot(rundata.lat,np.mean(psi[:,39,:],0))


plt.figure(2)
plt.contourf(rundata.lat,rundata.pfull,psi[210,:,:])
plt.show()
