#load in precip for a given run and produce a hovmoller plot at 15N

import sys
sys.path.append('/scratch/rg419/workdir_moist/python_bin')
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from time_av_xr import load_year_xr

def precip_hm_fn(inp_fol,nproc):

    #set run name
    run_fol = inp_fol + '/np' + nproc

    year = 1


    rundata = load_year_xr(run_fol, year)
    rundata.coords['pentad'] = (rundata.time // 5) + 1
    conv_p = rundata.convection_rain[:,37,:].groupby('pentad').mean(('time'))

    tot_p = xr.DataArray(np.zeros((72,128,5)), [('pentad', conv_p.pentad ), ('lon', rundata.lon), ('year', [2,3,4,5,6])])

    for year in range(2,7):
        rundata = load_year_xr(run_fol, year)
        rundata.coords['pentad'] = (rundata.time // 5) + 1
        #take time and zonal mean
        conv_p = rundata.convection_rain[:,37,:].groupby('pentad').mean(('time'))
        cond_p = rundata.condensation_rain[:,37,:].groupby('pentad').mean(('time'))
        tot_p[:,:,year-2] = (conv_p + cond_p)*86400

    tot_p_clim = tot_p.mean(('year'))

    #plot
    tot_p_clim.plot.contourf(x='lon', y='pentad',levels=np.arange(3.,31.,3.), add_label = False)
    plt.plot([float(rundata.lonb[21]),float(rundata.lonb[21])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[26]),float(rundata.lonb[26])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[30]),float(rundata.lonb[30])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[35]),float(rundata.lonb[35])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[39]),float(rundata.lonb[39])],[20.,60.],'k')
    plt.ylim(20,60)
    plt.xlim(60,180)
    plt.xlabel('Longitude')
    plt.ylabel('Pentad')
    plt.title('Precipitation at 15N, mm/day')
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/precip_hm.png')
    plt.clf()



def u_hm_fn(inp_fol,nproc):

    #set run name
    run_fol = inp_fol + '/np' + nproc
    year = 1
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['pentad'] = (rundata.time // 5) + 1
    conv_p = rundata.convection_rain[:,37,:].groupby('pentad').mean(('time'))

    u_low = xr.DataArray(np.zeros((72,128,5)), [('pentad', conv_p.pentad ), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    u_high = xr.DataArray(np.zeros((72,128,5)), [('pentad', conv_p.pentad ), ('lon', rundata.lon), ('year', [2,3,4,5,6])])

    for year in range(2,7):
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) + 1
        #take time and zonal mean
        u_low[:,:,year-2] = rundata.ucomp[:,2,37,:].groupby('pentad').mean(('time'))
        u_high[:,:,year-2] = rundata.ucomp[:,9,37,:].groupby('pentad').mean(('time'))


    u_low_clim = u_low.mean(('year'))
    u_high_clim = u_high.mean(('year'))

    #plot
    u_low_clim.plot.contourf(x='lon', y='pentad',levels=np.arange(-30.,31.,5.), add_label = False)
    plt.plot([float(rundata.lonb[21]),float(rundata.lonb[21])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[26]),float(rundata.lonb[26])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[30]),float(rundata.lonb[30])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[35]),float(rundata.lonb[35])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[39]),float(rundata.lonb[39])],[20.,60.],'k')
    plt.ylim(20,60)
    plt.xlim(60,180)
    plt.xlabel('Longitude')
    plt.ylabel('Pentad')
    plt.title('850 hPa zonal wind speed at 15N, m/s')
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/u850_hm.png')
    plt.clf()

    u_high_clim.plot.contourf(x='lon', y='pentad',levels=np.arange(-30.,31.,5.), add_label = False)
    plt.plot([float(rundata.lonb[21]),float(rundata.lonb[21])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[26]),float(rundata.lonb[26])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[30]),float(rundata.lonb[30])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[35]),float(rundata.lonb[35])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[39]),float(rundata.lonb[39])],[20.,60.],'k')
    plt.ylim(20,60)
    plt.xlim(60,180)
    plt.xlabel('Longitude')
    plt.ylabel('Pentad')
    plt.title('200 hPa zonal wind speed at 15N, m/s')
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/u200_hm.png')
    plt.clf()


def v_hm_fn(inp_fol,nproc):

    #set run name
    run_fol = inp_fol + '/np' + nproc
    year = 2
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['pentad'] = (rundata.time // 5) -71
    conv_p = rundata.vcomp.groupby('pentad').mean(('time'))

    v_low = xr.DataArray(np.zeros((72,128,5)), [('pentad', conv_p.pentad ), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    v_high = xr.DataArray(np.zeros((72,128,5)), [('pentad', conv_p.pentad ), ('lon', rundata.lon), ('year', [2,3,4,5,6])])

    for year in range(2,7):
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) + 1
        #take time and zonal mean
        v_low[:,:,year-2] = rundata.vcomp[:,2,37,:].groupby('pentad').mean(('time'))
        v_high[:,:,year-2] = rundata.vcomp[:,9,37,:].groupby('pentad').mean(('time'))


    v_low_clim = v_low.mean(('year'))
    v_high_clim = v_high.mean(('year'))

    #plot
    v_low_clim.plot.contourf(x='lon', y='pentad',levels=np.arange(-10.,11.,1.), add_label = False)
    plt.plot([float(rundata.lonb[21]),float(rundata.lonb[21])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[26]),float(rundata.lonb[26])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[30]),float(rundata.lonb[30])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[35]),float(rundata.lonb[35])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[39]),float(rundata.lonb[39])],[20.,60.],'k')
    plt.ylim(20,60)
    plt.xlim(60,180)
    plt.xlabel('Longitude')
    plt.ylabel('Pentad')
    plt.title('850 hPa meridional wind speed at 15N, m/s')
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/v850_hm.png')
    plt.clf()

    v_high_clim.plot.contourf(x='lon', y='pentad',levels=np.arange(-10.,11.,1.), add_label = False)
    plt.plot([float(rundata.lonb[21]),float(rundata.lonb[21])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[26]),float(rundata.lonb[26])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[30]),float(rundata.lonb[30])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[35]),float(rundata.lonb[35])],[20.,60.],'k')
    plt.plot([float(rundata.lonb[39]),float(rundata.lonb[39])],[20.,60.],'k')
    plt.ylim(20,60)
    plt.xlim(60,180)
    plt.xlabel('Longitude')
    plt.ylabel('Pentad')
    plt.title('200 hPa meridional wind speed at 15N, m/s')
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/v200_hm.png')
    plt.clf()


def ro_hm_fn(inp_fol,nproc):

    #set run name
    run_fol = inp_fol + '/np' + nproc
    year = 1
    rundata = load_year_xr(run_fol, year)
    rundata.coords['pentad'] = (rundata.time // 5) + 1
    conv_p = rundata.convection_rain[:,37,:].groupby('pentad').mean(('time'))

    vor = xr.DataArray(np.zeros((72,64,5)), [('pentad', conv_p.pentad ), ('lat', rundata.lat), ('year', [2,3,4,5,6])])

    for year in range(2,7):
        rundata = load_year_xr(run_fol, year)
        rundata.coords['pentad'] = (rundata.time // 5) + 1
        #take time and zonal mean
        vor[:,:,year-2] = rundata.vor[:,28,:,:].groupby('pentad').mean(('time','lon'))

    vor_clim = vor.mean(('year'))
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(rundata.lat *np.pi/180)
    ro_clim = -1.*vor_clim/f

    #plot
    ro_clim.plot.contourf(x='lat', y='pentad',levels=np.arange(0.,1.,0.1), add_label = False)
    plt.ylim(20,60)
    #plt.xlim(60,180)
    plt.xlabel('Latitude')
    plt.ylabel('Pentad')
    plt.title('206 hPa -vorticity/f')
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/ro206_hm.png')
    plt.clf()

