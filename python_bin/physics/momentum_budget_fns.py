#Load in u, v, omega, geopotential, evaluate the terms in the zonal momentum budget

from data_handling import load_year_xr
import numpy as np
import matplotlib.pyplot as plt
from finite_difference import cfd
import xarray as xr


def mombudg_2d_mn_fn(run_fol,nproc,years):
    year = years[0]
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['month'] = (rundata.time // 30) - 11
    mngrp = rundata.ucomp.groupby('month').mean(('time'))

    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(rundata.lat * np.pi/180)
    f = 2 * omega * np.sin(rundata.lat *np.pi/180)

    #Initialise arrays to load into
    u =     xr.DataArray(np.zeros((12,17,64,128,len(years))), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    v =     xr.DataArray(np.zeros((12,17,64,128,len(years))), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    w =     xr.DataArray(np.zeros((12,17,64,128,len(years))), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uu =    xr.DataArray(np.zeros((12,17,64,128,len(years))), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uv =    xr.DataArray(np.zeros((12,17,64,128,len(years))), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uw =    xr.DataArray(np.zeros((12,17,64,128,len(years))), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    phi =   xr.DataArray(np.zeros((12,17,64,128,len(years))), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    ddamp = xr.DataArray(np.zeros((12,17,64,128,len(years))), [('month', mngrp.month ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])

    for year in years:
        print year
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['month'] = (rundata.time // 30) + 1
        u[:,:,:,:,year-years[0]]     = rundata.ucomp.groupby('month').mean(('time'))
        v[:,:,:,:,year-years[0]]     = rundata.vcomp.groupby('month').mean(('time'))
        w[:,:,:,:,year-years[0]]     = rundata.omega.groupby('month').mean(('time'))
        uu[:,:,:,:,year-years[0]]    = rundata.ucomp_sq.groupby('month').mean(('time'))
        uv[:,:,:,:,year-years[0]]    = rundata.ucomp_vcomp.groupby('month').mean(('time'))
        uw[:,:,:,:,year-years[0]]    = rundata.ucomp_omega.groupby('month').mean(('time'))
        phi[:,:,:,:,year-years[0]]   = rundata.height.groupby('month').mean(('time'))
        ddamp[:,:,:,:,year-years[0]] = rundata.dt_ug_diffusion.groupby('month').mean(('time'))

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

    data_out = xr.Dataset({'dphidx_av': dphidx_av, 'fv_av': fv_av, 'ddamp_av': ddamp_av, 
                     'mom_eddy': mom_eddy, 'mom_mean': mom_mean, 'mom_sum': mom_sum, 'u_av': u_av, 'v_av': v_av  })

    return data_out





def mombudg_2d_pd_fn(run_fol,years, trange):
    year = years[0]
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['pentad'] = (rundata.time // 5) - 71
    mngrp = rundata.ucomp.groupby('pentad').mean(('time'))
    trl = trange[1]-trange[0]
    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(rundata.lat * np.pi/180)
    f = 2 * omega * np.sin(rundata.lat *np.pi/180)

    #Initialise arrays to load into
    u =     xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    v =     xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    w =     xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uu =    xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uv =    xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uw =    xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    phi =   xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    ddamp = xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])

    for year in years:
        print year
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) - 71
        u[:,:,:,:,year-years[0]]     = rundata.ucomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        v[:,:,:,:,year-years[0]]     = rundata.vcomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        w[:,:,:,:,year-years[0]]     = rundata.omega.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        uu[:,:,:,:,year-years[0]]    = rundata.ucomp_sq.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        uv[:,:,:,:,year-years[0]]    = rundata.ucomp_vcomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        uw[:,:,:,:,year-years[0]]    = rundata.ucomp_omega.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        phi[:,:,:,:,year-years[0]]   = rundata.height.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        ddamp[:,:,:,:,year-years[0]] = rundata.dt_ug_diffusion.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]

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
    dudx_av = xr.DataArray( cfd(u_av.values,u_av.lon*np.pi/180,3),   [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    ududx_av = u_av*dudx_av/coslat/a
    dudy_av = xr.DataArray( cfd((u_av*coslat).values,u_av.lat*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    vdudy_av = v_av*dudy_av/coslat/a
    dudp_av = xr.DataArray( cfd(u_av.values,u_av.pfull*100,1), [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    wdudp_av = w_av*dudp_av

    duueddx_av = xr.DataArray( cfd(uued_av.values,uued_av.lon*np.pi/180,3),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])
    duueddx_av = duueddx_av/coslat/a
    duveddy_av = xr.DataArray( cfd(uved_av.values ,uved_av.lat*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])
    duveddy_av = duveddy_av/coslat/coslat/a
    duweddp_av = xr.DataArray( cfd(uwed_av.values,uwed_av.pfull*100,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat), ('lon', uued_av.lon )])

    dphidx_av = xr.DataArray( cfd(phi_av.values,phi_av.lon*np.pi/180,3),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', phi_av.pfull ), ('lat', phi_av.lat), ('lon', phi_av.lon )])
    dphidx_av = -1*dphidx_av/coslat/a


    #evaluate sums of the terms
    mom_eddy = -(duueddx_av + duveddy_av + duweddp_av)
    mom_mean = -(ududx_av + vdudy_av + wdudp_av)
    mom_sum = fv_av + ddamp_av + dphidx_av + mom_eddy + mom_mean

    data_out = xr.Dataset({'dphidx_av': dphidx_av, 'fv_av': fv_av, 'ddamp_av': ddamp_av, 
                     'mom_eddy': mom_eddy, 'mom_mean': mom_mean, 'mom_sum': mom_sum , 'u_av': u_av, 'v_av': v_av   })


    return data_out



def mombudg_2d_an_fn(run_fol,years,era=True):

#inp_fol = experiment name
#nproc   = no of processors

    year = years[0]
    rundata = load_year_xr(run_fol, year, pinterp=True, era=era)


    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(rundata.lat * np.pi/180)
    f = 2 * omega * np.sin(rundata.lat *np.pi/180)
    if era==True:
        plevs = 17
    else:
        plevs=40
    #Initialise arrays to load into
    u =     xr.DataArray(np.zeros((plevs,64,128,len(years))), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    v =     xr.DataArray(np.zeros((plevs,64,128,len(years))), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    w =     xr.DataArray(np.zeros((plevs,64,128,len(years))), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uu =    xr.DataArray(np.zeros((plevs,64,128,len(years))), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uv =    xr.DataArray(np.zeros((plevs,64,128,len(years))), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uw =    xr.DataArray(np.zeros((plevs,64,128,len(years))), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    phi =   xr.DataArray(np.zeros((plevs,64,128,len(years))), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    ddamp = xr.DataArray(np.zeros((plevs,64,128,len(years))), [('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])

    for year in years:
        print year
        rundata = load_year_xr(run_fol, year, pinterp=True,era=era)
        rundata.coords['pentad'] = (rundata.time // 5) + 1
        u[:,:,:,year-years[0]]     = rundata.ucomp.mean(('time'))
        v[:,:,:,year-years[0]]     = rundata.vcomp.mean(('time'))
        w[:,:,:,year-years[0]]     = rundata.omega.mean(('time'))
        uu[:,:,:,year-years[0]]    = rundata.ucomp_sq.mean(('time'))
        uv[:,:,:,year-years[0]]    = rundata.ucomp_vcomp.mean(('time'))
        uw[:,:,:,year-years[0]]    = rundata.ucomp_omega.mean(('time'))
        phi[:,:,:,year-years[0]]   = rundata.height.mean(('time'))
        ddamp[:,:,:,year-years[0]] = rundata.dt_ug_diffusion.mean(('time'))

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


    data_out = xr.Dataset({'dphidx_av': dphidx_av, 'fv_av': fv_av, 'ddamp_av': ddamp_av, 
                     'mom_eddy': mom_eddy, 'mom_mean': mom_mean, 'mom_sum': mom_sum, 'u_av': u_av, 'v_av': v_av    })


    return data_out



   
if __name__ == "__main__":

        #set run name
        inp_fol = 'flat_10m_diags'
        nproc   = '16'
        years = [2,3,4,5,6]
	mombudg_2d_mn_fn(inp_fol,nproc,years)
