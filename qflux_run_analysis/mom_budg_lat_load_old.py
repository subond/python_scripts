# Evaluate and plot the momentum budget terms over a given latitude range

from data_handling import time_means, month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from finite_difference import cfd
import gc

plt.rc('text', usetex=True)
font = {'size'   : 18}
plt.rc('font', **font)
mathtxt={'fontset':'custom', 'default':'regular'}
plt.rc('mathtext',**mathtxt)

def load_mom_vars(run, output_file):
    data = time_means(run, [121,481], filename='plev_pentad', timeav='pentad')
    print 'data passed'

    mom_budg = xr.Dataset({'ucomp': data.ucomp, 'vcomp': data.vcomp, 'omega': data.omega, 
                     'ucomp_sq': data.ucomp_sq, 'ucomp_vcomp': data.ucomp_vcomp, 'ucomp_omega': data.ucomp_omega, 
                     'height': data.height, 'dt_ug_diffusion': data.dt_ug_diffusion  })
    print 'mom_budg defined'

    mom_budg.to_netcdf(path=output_file, mode='w')
    print 'mom_budg saved'

    data.close()
    mom_budg.close()
    print 'data closed'


load_mom_vars('q_flux_test','q_flux_test.nc')
print 'qflux done'