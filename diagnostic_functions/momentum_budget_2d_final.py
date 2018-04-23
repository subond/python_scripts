#Load in u, v, omega, geopotential, evaluate the terms in the zonal momentum budget

import sys
sys.path.append('/scratch/rg419/workdir_moist/python_bin')
from time_av_xr import load_year_xr
import numpy as np
import matplotlib.pyplot as plt
from finite_difference import cfd
import xarray as xr


def mombudg_2d_mn_fn(inp_fol,nproc):

#inp_fol = experiment name
#nproc   = no of processors

    #set run name
    run_fol = inp_fol + '/np' + nproc
    year = 2
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['month'] = (rundata.time // 30) - 11
    mngrp = rundata.ucomp.groupby('month').mean(('time'))

    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(rundata.lat * np.pi/180)
    f = 2 * omega * np.sin(rundata.lat *np.pi/180)

    #Initialise arrays to load into
    u =     xr.DataArray(np.zeros((12,17,64,128,5)), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    v =     xr.DataArray(np.zeros((12,17,64,128,5)), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    w =     xr.DataArray(np.zeros((12,17,64,128,5)), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    uu =    xr.DataArray(np.zeros((12,17,64,128,5)), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    uv =    xr.DataArray(np.zeros((12,17,64,128,5)), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    uw =    xr.DataArray(np.zeros((12,17,64,128,5)), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    phi =   xr.DataArray(np.zeros((12,17,64,128,5)), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    ddamp = xr.DataArray(np.zeros((12,17,64,128,5)), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])

    for year in range(2,7):
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['month'] = (rundata.time // 30) + 1
        u[:,:,:,:,year-2]     = rundata.ucomp.groupby('month').mean(('time'))
        v[:,:,:,:,year-2]     = rundata.vcomp.groupby('month').mean(('time'))
        w[:,:,:,:,year-2]     = rundata.omega.groupby('month').mean(('time'))
        uu[:,:,:,:,year-2]    = rundata.ucomp_sq.groupby('month').mean(('time'))
        uv[:,:,:,:,year-2]    = rundata.ucomp_vcomp.groupby('month').mean(('time'))
        uw[:,:,:,:,year-2]    = rundata.ucomp_omega.groupby('month').mean(('time'))
        phi[:,:,:,:,year-2]   = rundata.height.groupby('month').mean(('time'))
        ddamp[:,:,:,:,year-2] = rundata.dt_ug_diffusion.groupby('month').mean(('time'))

    #Take averages: u, v, uv, w, phi, damping terms
    u_av = u.mean(('year'))
    v_av = v.mean(('year'))
    w_av = w.mean(('year'))
    uu_av = uu.mean(('year'))
    uv_av = uv.mean(('year'))
    uw_av = uw.mean(('year'))
    phi_av = phi.mean(('year'))*9.8

    ddamp_av = ddamp.mean(('year'))
    fv_av = v_av*f 

    #Evaluate mean eddy products
    uued_av  =  uu_av - u_av*u_av
    uved_av  = (uv_av - u_av*v_av)*coslat*coslat
    uwed_av  =  uw_av - u_av*w_av

    #Evaluate gradients needed
    dudx_av = xr.DataArray( cfd(u_av.values,u_av.lon*np.pi/180,3),   [('month', mngrp.month ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    ududx_av = u_av*dudx_av/coslat/a
    dudy_av = xr.DataArray( cfd((u_av*coslat).values,u_av.lat*np.pi/180,2),   [('month', mngrp.month ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    vdudy_av = v_av*dudy_av/coslat/a
    dudp_av = xr.DataArray( cfd(u_av.values,u_av.pfull*100,1), [('month', mngrp.month ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    wdudp_av = w_av*dudp_av

    duueddx_av = xr.DataArray( cfd(uued_av.values,uued_av.lon*np.pi/180,3),   [('month', mngrp.month ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])
    duueddx_av = duueddx_av/coslat/a
    duveddy_av = xr.DataArray( cfd(uved_av.values ,uved_av.lat*np.pi/180,2),   [('month', mngrp.month ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])
    duveddy_av = duveddy_av/coslat/coslat/a
    duweddp_av = xr.DataArray( cfd(uwed_av.values,uwed_av.pfull*100,1),   [('month', mngrp.month ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])

    dphidx_av = xr.DataArray( cfd(phi_av.values,phi_av.lon*np.pi/180,3),   [('month', mngrp.month ), ('pfull', phi_av.pfull ), ('lat', phi_av.lat), ('lon', phi_av.lon )])
    dphidx_av = -1*dphidx_av/coslat/a


    #evaluate sums of the terms
    mom_eddy = -(duueddx_av + duveddy_av + duweddp_av)
    mom_mean = -(ududx_av + vdudy_av + wdudp_av)
    mom_sum = fv_av + ddamp_av + dphidx_av + mom_eddy + mom_mean

    #Plot up seasonal cycle for each term at 200 and 850hPa and save

    dphidx_av[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/dphidx_200.png')
    fv_av[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/fv_200.png')
    mom_eddy[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/momeddy_200.png')
    mom_mean[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/mommean_200.png')
    ddamp_av[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/ddamp_200.png')
    mom_sum[:,9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/momsum_200.png')

    dphidx_av[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/dphidx_850.png')
    fv_av[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/fv_850.png')
    mom_eddy[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/momeddy_850.png')
    mom_mean[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/mommean_850.png')
    ddamp_av[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/ddamp_850.png')
    mom_sum[:,2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001), col='month', col_wrap=4)
    plt.ylim(-30,30)
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/momsum_850.png')

    data_out = xr.Dataset({'dphidx_av': dphidx_av, 'fv_av': fv_av, 'ddamp_av': ddamp_av, 
                     'mom_eddy': mom_eddy, 'mom_mean': mom_mean, 'mom_sum': mom_sum  })

    return data_out




def mombudg_2d_pd_fn(inp_fol,nproc):

#inp_fol = experiment name
#nproc   = no of processors

    #set run name
    run_fol = inp_fol + '/np' + nproc
    year = 2
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['pentad'] = (rundata.time // 5) - 71
    mngrp = rundata.ucomp.groupby('pentad').mean(('time'))

    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(rundata.lat * np.pi/180)
    f = 2 * omega * np.sin(rundata.lat *np.pi/180)

    #Initialise arrays to load into
    u =     xr.DataArray(np.zeros((72,17,64,128,5)), [('pentad', mngrp.pentad ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    v =     xr.DataArray(np.zeros((72,17,64,128,5)), [('pentad', mngrp.pentad ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    w =     xr.DataArray(np.zeros((72,17,64,128,5)), [('pentad', mngrp.pentad ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    uu =    xr.DataArray(np.zeros((72,17,64,128,5)), [('pentad', mngrp.pentad ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    uv =    xr.DataArray(np.zeros((72,17,64,128,5)), [('pentad', mngrp.pentad ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    uw =    xr.DataArray(np.zeros((72,17,64,128,5)), [('pentad', mngrp.pentad ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    phi =   xr.DataArray(np.zeros((72,17,64,128,5)), [('pentad', mngrp.pentad ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    ddamp = xr.DataArray(np.zeros((72,17,64,128,5)), [('pentad', mngrp.pentad ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])

    for year in range(2,7):
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) - 71
        u[:,:,:,:,year-2]     = rundata.ucomp.groupby('pentad').mean(('time'))
        v[:,:,:,:,year-2]     = rundata.vcomp.groupby('pentad').mean(('time'))
        w[:,:,:,:,year-2]     = rundata.omega.groupby('pentad').mean(('time'))
        uu[:,:,:,:,year-2]    = rundata.ucomp_sq.groupby('pentad').mean(('time'))
        uv[:,:,:,:,year-2]    = rundata.ucomp_vcomp.groupby('pentad').mean(('time'))
        uw[:,:,:,:,year-2]    = rundata.ucomp_omega.groupby('pentad').mean(('time'))
        phi[:,:,:,:,year-2]   = rundata.height.groupby('pentad').mean(('time'))
        ddamp[:,:,:,:,year-2] = rundata.dt_ug_diffusion.groupby('pentad').mean(('time'))

    #Take averages: u, v, uv, w, phi, damping terms
    u_av = u.mean(('year'))
    v_av = v.mean(('year'))
    w_av = w.mean(('year'))
    uu_av = uu.mean(('year'))
    uv_av = uv.mean(('year'))
    uw_av = uw.mean(('year'))
    phi_av = phi.mean(('year'))*9.8

    ddamp_av = ddamp.mean(('year'))
    fv_av = v_av*f 

    #Evaluate mean eddy products
    uued_av  =  uu_av - u_av*u_av
    uved_av  = (uv_av - u_av*v_av)*coslat*coslat
    uwed_av  =  uw_av - u_av*w_av

    #Evaluate gradients needed
    dudx_av = xr.DataArray( cfd(u_av.values,u_av.lon*np.pi/180,3),   [('pentad', mngrp.pentad ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    ududx_av = u_av*dudx_av/coslat/a
    dudy_av = xr.DataArray( cfd((u_av*coslat).values,u_av.lat*np.pi/180,2),   [('pentad', mngrp.pentad ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    vdudy_av = v_av*dudy_av/coslat/a
    dudp_av = xr.DataArray( cfd(u_av.values,u_av.pfull*100,1), [('pentad', mngrp.pentad ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    wdudp_av = w_av*dudp_av

    duueddx_av = xr.DataArray( cfd(uued_av.values,uued_av.lon*np.pi/180,3),   [('pentad', mngrp.pentad ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])
    duueddx_av = duueddx_av/coslat/a
    duveddy_av = xr.DataArray( cfd(uved_av.values ,uved_av.lat*np.pi/180,2),   [('pentad', mngrp.pentad ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])
    duveddy_av = duveddy_av/coslat/coslat/a
    duweddp_av = xr.DataArray( cfd(uwed_av.values,uwed_av.pfull*100,1),   [('pentad', mngrp.pentad ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])

    dphidx_av = xr.DataArray( cfd(phi_av.values,phi_av.lon*np.pi/180,3),   [('pentad', mngrp.pentad ), ('pfull', phi_av.pfull ), ('lat', phi_av.lat), ('lon', phi_av.lon )])
    dphidx_av = -1*dphidx_av/coslat/a


    #evaluate sums of the terms
    mom_eddy = -(duueddx_av + duveddy_av + duweddp_av)
    mom_mean = -(ududx_av + vdudy_av + wdudp_av)
    mom_sum = fv_av + ddamp_av + dphidx_av + mom_eddy + mom_mean

    data_out = xr.Dataset({'dphidx_av': dphidx_av, 'fv_av': fv_av, 'ddamp_av': ddamp_av, 
                     'mom_eddy': mom_eddy, 'mom_mean': mom_mean, 'mom_sum': mom_sum  })


    return data_out



def mombudg_2d_an_fn(inp_fol,nproc):

#inp_fol = experiment name
#nproc   = no of processors

    #set run name
    run_fol = inp_fol + '/np' + nproc
    year = 2
    rundata = load_year_xr(run_fol, year, pinterp=True)


    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(rundata.lat * np.pi/180)
    f = 2 * omega * np.sin(rundata.lat *np.pi/180)

    #Initialise arrays to load into
    u =     xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    v =     xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    w =     xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    uu =    xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    uv =    xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    uw =    xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    phi =   xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])
    ddamp = xr.DataArray(np.zeros((17,64,128,5)), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', [2,3,4,5,6])])

    for year in range(2,7):
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) + 1
        u[:,:,:,year-2]     = rundata.ucomp.mean(('time'))
        v[:,:,:,year-2]     = rundata.vcomp.mean(('time'))
        w[:,:,:,year-2]     = rundata.omega.mean(('time'))
        uu[:,:,:,year-2]    = rundata.ucomp_sq.mean(('time'))
        uv[:,:,:,year-2]    = rundata.ucomp_vcomp.mean(('time'))
        uw[:,:,:,year-2]    = rundata.ucomp_omega.mean(('time'))
        phi[:,:,:,year-2]   = rundata.height.mean(('time'))
        ddamp[:,:,:,year-2] = rundata.dt_ug_diffusion.mean(('time'))

    #Take averages: u, v, uv, w, phi, damping terms
    u_av = u.mean(('year'))
    v_av = v.mean(('year'))
    w_av = w.mean(('year'))
    uu_av = uu.mean(('year'))
    uv_av = uv.mean(('year'))
    uw_av = uw.mean(('year'))
    phi_av = phi.mean(('year'))*9.8

    ddamp_av = ddamp.mean(('year'))
    fv_av = v_av*f 

    #Evaluate mean eddy products
    uued_av  =  uu_av - u_av*u_av
    uved_av  = (uv_av - u_av*v_av)*coslat*coslat
    uwed_av  =  uw_av - u_av*w_av

    #Evaluate gradients needed
    dudx_av = xr.DataArray( cfd(u_av.values,u_av.lon*np.pi/180,2),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    ududx_av = u_av*dudx_av/coslat/a
    dudy_av = xr.DataArray( cfd((u_av*coslat).values,u_av.lat*np.pi/180,1),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    vdudy_av = v_av*dudy_av/coslat/a
    dudp_av = xr.DataArray( cfd(u_av.values,u_av.pfull*100,0), [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    wdudp_av = w_av*dudp_av

    duueddx_av = xr.DataArray( cfd(uued_av.values,uued_av.lon*np.pi/180,2),   [('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])
    duueddx_av = duueddx_av/coslat/a
    duveddy_av = xr.DataArray( cfd(uved_av.values ,uved_av.lat*np.pi/180,1),   [('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])
    duveddy_av = duveddy_av/coslat/coslat/a
    duweddp_av = xr.DataArray( cfd(uwed_av.values,uwed_av.pfull*100,0),   [('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])

    dphidx_av = xr.DataArray( cfd(phi_av.values,phi_av.lon*np.pi/180,2),   [('pfull', phi_av.pfull ), ('lat', phi_av.lat), ('lon', phi_av.lon )])
    dphidx_av = -1*dphidx_av/coslat/a


    #evaluate sums of the terms
    mom_eddy = -(duueddx_av + duveddy_av + duweddp_av)
    mom_mean = -(ududx_av + vdudy_av + wdudp_av)
    mom_sum = fv_av + ddamp_av + dphidx_av + mom_eddy + mom_mean


    #Plot up seasonal cycle for each term at 200 and 850hPa and save

    dphidx_av[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_dphidx_200.png')
    plt.clf()
    fv_av[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_fv_200.png')
    plt.clf()
    mom_eddy[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_momeddy_200.png')
    plt.clf()
    mom_mean[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_mommean_200.png')
    plt.clf()
    ddamp_av[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_ddamp_200.png')
    plt.clf()
    mom_sum[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.00021,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_momsum_200.png')
    plt.clf()

    dphidx_av[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_dphidx_850.png')
    plt.clf()
    fv_av[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_fv_850.png')
    plt.clf()
    mom_eddy[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_momeddy_850.png')
    plt.clf()
    mom_mean[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_mommean_850.png')
    plt.clf()
    ddamp_av[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_ddamp_850.png')
    plt.clf()
    mom_sum[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.00011,0.00001))
    plt.savefig('/scratch/rg419/workdir_moist/'+inp_fol+'/an_momsum_850.png')
    plt.clf()



    data_out = xr.Dataset({'dphidx_av': dphidx_av, 'fv_av': fv_av, 'ddamp_av': ddamp_av, 
                     'mom_eddy': mom_eddy, 'mom_mean': mom_mean, 'mom_sum': mom_sum  })


    return data_out



   
if __name__ == "__main__":

        #set run name
        inp_fol = 'flat_10m_diags'
        nproc   = '16'

	mombudg_2d_an_fn(inp_fol,nproc)
	mombudg_2d_mn_fn(inp_fol,nproc)
	#mombudg_2d_pd_fn(inp_fol,nproc)
