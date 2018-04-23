# Produce B+S style plots of zonal mean precip over time for the two aquaplanet runs.

from plotting import load_bs_vars
from data_handling import pentad_dic, month_dic
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

#plt.rc('text', usetex=True)
font = {#'family' : 'sans-serif','sans-serif':['Helvetica'],
        #'weight' : 'bold',
        'size'   : 18}

plt.rc('font', **font)


def bs_plot(inp_fol, years):
    data = load_bs_vars(inp_fol, years)
    #pd_dic = pentad_dic(1)
    mn_dic = month_dic(1)
    #tickspace = range(9,80,18)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    ax=data.totp.plot.contourf(x='pentad', y='lat',levels=np.arange(6.,31.,3.), add_label = False, add_colorbar=False)
    cb1=plt.colorbar(ax)
    cb1.set_label('Precip, mm/day')
    cs = data.t_surf.plot.contour(x='pentad', y='lat',levels=range(250,321,10), add_label = False, colors='w', add_colorbar=False)
    plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
    plt.ylim((-45,45))
    plt.xlim((1,73))
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.ylabel('Latitude')
    #plt.title('Zonal mean precipitation evolution, mm/day')
    plt.tight_layout()
    
    
bs_plot('aquaplanet_10m', range(11,41))
plt.savefig('/scratch/rg419/plots/mom_budg_work/precip_t_10m.png')
plt.clf()
bs_plot('aquaplanet_2m', range(11,41))
plt.savefig('/scratch/rg419/plots/mom_budg_work/precip_t_2m.png')
plt.clf()