# Plot the terms in the momentum budget for the topography, flat, and aquamountain runs as lon-time plots based around 15N across the monsoon region. May want to play with the range of latitudes averaged across.

from physics import mombudg_lon_pd_fn
from data_handling import load_year_xr
from pentads import pentad_dic, month_dic
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


plt.rc('text', usetex=True)
font = {'family' : 'sans-serif','sans-serif':['Helvetica'],
        'weight' : 'bold',
        'size'   : 18}

plt.rc('font', **font)

def plot_pd_mom_fn(inp_fol, years):
    
    year = years[0]
    rundata = load_year_xr(inp_fol, year)
    tot_p = xr.DataArray(np.zeros((72,128,len(years))), [('pentad', range(1,73) ), ('lon', rundata.lon), ('year', years)])

    i=0
    for year in years:
        print year
        rundata = load_year_xr(inp_fol, year)
        rundata.coords['pentad'] = (rundata.time // 5) -71
        conv_p = rundata.convection_rain[:,37,:].groupby('pentad').mean(('time'))
        cond_p = rundata.condensation_rain[:,37,:].groupby('pentad').mean(('time'))
        tot_p[:,:,i] = (conv_p + cond_p)*86400
        i=i+1

    tot_p_clim = tot_p.mean(('year'))
    
    data,wspd = mombudg_lon_pd_fn(inp_fol, years, [19,61])

    mn_dic = month_dic(1)
    tickspace = range(25,60,6)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    def plot_mom_var(var,levels):
        var_dic = {'fv_av': 'fv',
                    'mom_eddy': 'Eddy advective terms',
                    'mom_mean': 'Mean state advective terms',
                    'mom_sum': 'Residual',
                    'dphidx_av': 'Geopotential gradient'}
        if var == 'fv_mn_imb':
            plot_data = data.data_vars['fv_av'] + data.data_vars['mom_mean']
        else:
            plot_data = data.data_vars[var]
        plot_data = plot_data*10000.
        g=plot_data.plot.contourf(x='lon', y='pentad',levels=levels, add_labels = False, add_colorbar=False,extend='both')
        cb1=plt.colorbar(g)
        cb1.set_label('$\displaystyle10^{-4}m/s^2$')
        plt.quiver(data.lon[::3], data.pentad[::2], wspd.u_av[::2,::3], wspd.v_av[::2,::3], headwidth=3)#, scale=500.)
        #tickspace = range(25,60,10)
        #labels = [pd_dic[k] for k in tickspace]
        plt.ylim(20,60)
        plt.xlim(60,180)
        if not 'aqua' in inp_fol:
            plt.plot([float(rundata.lonb[21]),float(rundata.lonb[21])],[0.,60.],'k')
            plt.plot([float(rundata.lonb[26]),float(rundata.lonb[26])],[0.,60.],'k')
            plt.plot([float(rundata.lonb[30]),float(rundata.lonb[30])],[0.,60.],'k')
            plt.plot([float(rundata.lonb[35]),float(rundata.lonb[35])],[0.,60.],'k')
            plt.plot([float(rundata.lonb[39]),float(rundata.lonb[39])],[0.,60.],'k')
            plt.fill_between([float(rundata.lonb[26]),float(rundata.lonb[30])], 0., 60., facecolor='gray', alpha=0.5)
            plt.fill_between([float(rundata.lonb[35]),float(rundata.lonb[39])], 0., 60., facecolor='gray', alpha=0.5)
        else:
            plt.plot([float(rundata.lonb[21]),float(rundata.lonb[21])],[0.,60.],'k--')
            plt.plot([float(rundata.lonb[26]),float(rundata.lonb[26])],[0.,60.],'k--')
            plt.plot([float(rundata.lonb[30]),float(rundata.lonb[30])],[0.,60.],'k--')
            plt.plot([float(rundata.lonb[35]),float(rundata.lonb[35])],[0.,60.],'k--')
            plt.plot([float(rundata.lonb[39]),float(rundata.lonb[39])],[0.,60.],'k--')
        tot_p_clim.plot.contour(x='lon', y='pentad', colors='k', levels=[-1000.,6.,1000.], add_labels = False, add_colorbar=False)
        #plt.contour(tot_p_clim.lon, tot_p_clim.pentad, tot_p_clim.values, colors='k', levels=[-1000.,6.,1000.], linewidth=10.)
        plt.xlabel('Longitude')
        plt.yticks(tickspace,labels,rotation=25)
        #plt.ylabel('Pentad')
        #plt.title(var_dic[var] + ', $\displaystyle10^{-4}m/s^2$')
        plt.tight_layout()
        plt.savefig('/scratch/rg419/plots/mom_budg_work/'+inp_fol+'/'+var+'_30yr.png')
        plt.close()
        
    vars = ['fv_mn_imb','fv_av','mom_eddy','mom_mean','mom_sum','dphidx_av']
    for var in vars:
        plot_mom_var(var,np.arange(-2.,2.1,0.2))            
        #plot_mom_var(var,9,np.arange(-0.0002,0.00021,0.00002))            

plot_pd_mom_fn('aquamountain_10m',range(11,41))
#plot_pd_mom_fn('aquaplanet_2m',range(11,41))
plot_pd_mom_fn('flat_10m',range(11,41))
plot_pd_mom_fn('topo_10m',range(11,41))

