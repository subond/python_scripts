"""
Take divergence of u,v,w and plot hm - are we incompressible?? (01/09/2018)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, model_constants as mc
from pylab import rcParams

plot_dir = '/scratch/rg419/plots/paper_2_figs/revisions/'

def incompressibility(run, lev=850.):

    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    rho = data.pfull*100./data.temp/mc.rdgas
    rho = data.pfull*100./mc.rdgas/data.temp / (1 + (mc.rvgas - mc.rdgas)/mc.rdgas*data.sphum)
    
    dudx = gr.ddx(data.ucomp * rho)
    dvdy = gr.ddy(data.vcomp * rho)
    dwdp = gr.ddp(data.omega * rho)
    
    divU = (dudx + dvdy + dwdp).mean('lon').sel(pfull=lev) * 86400.
    
    rcParams['figure.figsize'] = 5.5, 3.
    rcParams['font.size'] = 14
    
    fig, ax = plt.subplots(1, sharex=True)
    
    divU.plot.contourf(ax=ax, x='xofyear', y='lat', add_labels=False)
    
    ax.set_xlabel('Pentad')
    ax.set_ylabel('Latitude')
    
    plt.subplots_adjust(left=0.15, right=0.95, top=0.97, bottom=0.2)
    
    plt.savefig(plot_dir+'incompressibility_rho_' + run + '.pdf', format='pdf')
    plt.close()
    
incompressibility('sn_1.000_evap_fluxes_heattrans')