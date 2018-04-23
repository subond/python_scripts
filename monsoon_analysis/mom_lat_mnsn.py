# Make contour plots of the terms in the momentum budget at 200hPa plus the wind vectors for the topography run over the onset period. 

from physics import mombudg_lev_pd_fn
from data_handling import pentad_dic, month_dic
from plotting import load_bs_vars
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


plt.rc('text', usetex=True)
font = {'size'   : 18}
plt.rc('font', **font)
mathtxt={'fontset':'custom', 'default':'regular'}
plt.rc('mathtext',**mathtxt)


def plot_pd_mom_fn(inp_fol):
    
    data_p = load_bs_vars(inp_fol, range(11,41))
    data_file = '/scratch/rg419/plots/monsoon_analysis/'+inp_fol+'/mombudg_data.nc'
    data= xr.open_dataset( data_file) 
    print 'data loaded, starting plotting'
    land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/topo_10m/input/land.nc'
    land = xr.open_dataset( land_file)
    land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
    mn_dic = month_dic(1)
    tickspace = range(13,72,18)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
        
    def plot_mom_var(var,levels):

        if var == 'fv_mn_imb':
            plot_dataa = data.data_vars['fv_av'] + data.data_vars['mom_mean']
        else:
            plot_dataa = data.data_vars[var]
        plot_data = plot_dataa[:,:,22:65].mean(('lon'))*10000.
        
        g=plot_data.plot.contourf(x='pentad', y='lat',levels=levels, add_label = False, add_colorbar=False, extend='both')
        cb1=plt.colorbar(g)
        cb1.set_label('$\displaystyle10^{-4}m/s^2$')
        cs=data_p.totp.plot.contour(x='pentad', y='lat',levels=np.arange(6.,31.,6.), add_label = False, add_colorbar=False,colors='k')
        plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
        plt.ylim(0,90)
        plt.xlim(1,73)
        plt.xlabel('')
        plt.xticks(tickspace,labels,rotation=25)
        plt.ylabel('Latitude')
        #plt.plot([float(data.pentad[36]),float(data.pentad[36])],[0.,90.],'k--')
        plt.tight_layout()
        plt.title('')
        #plt.title(var_dic[var] + ', $\displaystyle10^{-4}m/s^2$')
        plt.savefig('/scratch/rg419/plots/monsoon_analysis/'+inp_fol+'/'+var+'_lattime.png')
        plt.close()
    
        
        
    vars = ['dphidx_av','fv_av','mom_eddy','mom_mean','mom_sum','fv_mn_imb']
    for var in vars:
        plot_mom_var(var,np.arange(-2.,2.1,0.4))            

plot_pd_mom_fn('aquaplanet_10m')
print 'aquaplanet_10m done'
plot_pd_mom_fn('aquaplanet_2m')
print 'aquaplanet_2m done'
plot_pd_mom_fn('flat_10m')
print 'flat_10m done'
plot_pd_mom_fn('aquamountain_10m')
print 'aquamountain_10m done'
plot_pd_mom_fn('topo_10m')
print 'topo_10m done'
