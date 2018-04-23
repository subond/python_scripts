# Produce plots of precip and wind vectors at 15N over time for the monsoon region in the topography run for presentation/paper
# Also produce same for flat and aquamountain runs for comparison.

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from data_handling import load_year_xr, pentad_dic

def precip_hm(run_fol,years):

    year = years[0]
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['pentad'] = (rundata.time // 5) -71
    mngrp = rundata.ucomp.groupby('pentad').mean(('time'))

    tot_p = xr.DataArray(np.zeros((72,128,len(years))), [('pentad', range(1,73) ), ('lon', rundata.lon), ('year', years)])
    u = xr.DataArray(np.zeros((72,17,128,len(years))), [('pentad', range(1,73)  ), ('pfull', rundata.pfull), ('lon', rundata.lon), ('year', years)])
    v = xr.DataArray(np.zeros((72,17,128,len(years))), [('pentad', range(1,73)  ), ('pfull', rundata.pfull), ('lon', rundata.lon), ('year', years)])
    phi = xr.DataArray(np.zeros((72,17,128,len(years))), [('pentad', range(1,73)  ), ('pfull', rundata.pfull), ('lon', rundata.lon), ('year', years)])

    i=0
    for year in years:
        print year
        rundata = load_year_xr(run_fol, year)
        rundata.coords['pentad'] = (rundata.time // 5) -71
        #take time mean
        conv_p = rundata.convection_rain[:,37,:].groupby('pentad').mean(('time'))
        cond_p = rundata.condensation_rain[:,37,:].groupby('pentad').mean(('time'))
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) -71
        u[:,:,:,i] = rundata.ucomp[:,:,37,:].groupby('pentad').mean(('time'))
        v[:,:,:,i] = rundata.vcomp[:,:,37,:].groupby('pentad').mean(('time'))
        phi[:,:,:,i] = rundata.height[:,:,37,:].groupby('pentad').mean(('time'))
        tot_p[:,:,i] = (conv_p + cond_p)*86400
        i=i+1

    tot_p_clim = tot_p.mean(('year'))
    u_clim = u.mean(('year'))
    v_clim = v.mean(('year'))
    phi_clim = phi.mean(('year'))
    phi_ed_clim = phi_clim - phi_clim.mean(('lon'))

    pd_dic = pentad_dic(1)
    tickspace = range(20,61,5)
    labels = [pd_dic[k] for k in tickspace]  
    
    #plot    
    plt.figure(0)
    phi_ed_clim[:,2,:].plot.contourf(x='lon', y='pentad', levels=np.arange(-40.,40.,2.), add_label = False)
    tot_p_clim.plot.contour(x='lon', y='pentad',levels=[-1000.,6.,1000.], add_label = False, add_colorbar=False, colors='k')
    plt.yticks(tickspace,labels,rotation=25)
    plt.quiver(u_clim.lon[::3], u_clim.pentad[::2], u_clim[::2,2,::3], v_clim[::2,2,::3], headwidth=3)#, scale=500.)
    plt.plot([float(rundata.lonb[21]),float(rundata.lonb[21])],[0.,60.],'k')
    plt.plot([float(rundata.lonb[26]),float(rundata.lonb[26])],[0.,60.],'k')
    plt.plot([float(rundata.lonb[30]),float(rundata.lonb[30])],[0.,60.],'k')
    plt.plot([float(rundata.lonb[35]),float(rundata.lonb[35])],[0.,60.],'k')
    plt.plot([float(rundata.lonb[39]),float(rundata.lonb[39])],[0.,60.],'k')
    plt.ylim(20,60)
    plt.xlim(60,180)
    plt.xlabel('Longitude')
    plt.ylabel('Pentad')
    plt.title('850hPa geopotential anomalies, m')
    plt.savefig('/scratch/rg419/plots/mom_budg_work/geopot_hm_850.png')
    plt.clf()

    plt.figure(0)
    phi_ed_clim[:,9,:].plot.contourf(x='lon', y='pentad', levels=np.arange(-40.,40.,2.), add_label = False)
    tot_p_clim.plot.contour(x='lon', y='pentad',levels=[-1000.,6.,1000.], add_label = False, add_colorbar=False, colors='k')
    plt.yticks(tickspace,labels,rotation=25)
    plt.quiver(u_clim.lon[::3], u_clim.pentad[::2], u_clim[::2,9,::3], v_clim[::2,9,::3], headwidth=3)#, scale=500.)
    plt.plot([float(rundata.lonb[21]),float(rundata.lonb[21])],[0.,60.],'k')
    plt.plot([float(rundata.lonb[26]),float(rundata.lonb[26])],[0.,60.],'k')
    plt.plot([float(rundata.lonb[30]),float(rundata.lonb[30])],[0.,60.],'k')
    plt.plot([float(rundata.lonb[35]),float(rundata.lonb[35])],[0.,60.],'k')
    plt.plot([float(rundata.lonb[39]),float(rundata.lonb[39])],[0.,60.],'k')
    plt.ylim(20,60)
    plt.xlim(60,180)
    plt.xlabel('Longitude')
    plt.ylabel('Pentad')
    plt.title('200hPa geopotential anomalies, m')
    plt.savefig('/scratch/rg419/plots/mom_budg_work/geopot_hm_200.png')
    plt.clf()
    
if __name__ == "__main__":

    precip_hm('topo_10m',range(21,41))

