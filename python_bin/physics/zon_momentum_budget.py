#Load in u, v, omega, geopotential, evaluate the terms in the zonal momentum budget

from data_handling import load_year_xr
import numpy as np
import matplotlib.pyplot as plt
from finite_difference import cfd
import xarray as xr


def mombudg_pd_fn(run_fol,years):
#inp_fol = experiment name
#nproc   = no of processors

    year = years[0]
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['pentad'] = (rundata.time // 5) - 71

    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(rundata.lat * np.pi/180)
    f = 2 * omega * np.sin(rundata.lat *np.pi/180)

    #Initialise arrays to load into
    u =     xr.DataArray(np.zeros((72,17,64,len(years))), [('pentad', range(1,73) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('year', years )])
    v =     xr.DataArray(np.zeros((72,17,64,len(years))), [('pentad', range(1,73) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('year', years )])
    w =     xr.DataArray(np.zeros((72,17,64,len(years))), [('pentad', range(1,73) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('year', years )])
    uu =    xr.DataArray(np.zeros((72,17,64,len(years))), [('pentad', range(1,73) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('year', years )])
    uv =    xr.DataArray(np.zeros((72,17,64,len(years))), [('pentad', range(1,73) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('year', years )])
    uw =    xr.DataArray(np.zeros((72,17,64,len(years))), [('pentad', range(1,73) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('year', years )])
    phi =   xr.DataArray(np.zeros((72,17,64,len(years))), [('pentad', range(1,73) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('year', years )])
    ddamp = xr.DataArray(np.zeros((72,17,64,len(years))), [('pentad', range(1,73) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('year', years )])

    for year in years:
        print year
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) - 71
        u[:,:,:,year-years[0]]     = rundata.ucomp.groupby('pentad').mean(('time','lon'))
        v[:,:,:,year-years[0]]     = rundata.vcomp.groupby('pentad').mean(('time','lon'))
        w[:,:,:,year-years[0]]     = rundata.omega.groupby('pentad').mean(('time','lon'))
        uu[:,:,:,year-years[0]]    = rundata.ucomp_sq.groupby('pentad').mean(('time','lon'))
        uv[:,:,:,year-years[0]]    = rundata.ucomp_vcomp.groupby('pentad').mean(('time','lon'))
        uw[:,:,:,year-years[0]]    = rundata.ucomp_omega.groupby('pentad').mean(('time','lon'))
        phi[:,:,:,year-years[0]]   = rundata.height.groupby('pentad').mean(('time','lon'))
        ddamp[:,:,:,year-years[0]] = rundata.dt_ug_diffusion.groupby('pentad').mean(('time','lon'))

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

    dudy_av = xr.DataArray( cfd((u_av*coslat).values,u_av.lat*np.pi/180,2),   [('pentad', range(1,73)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
    vdudy_av = v_av*dudy_av/coslat/a
    dudp_av = xr.DataArray( cfd(u_av.values,u_av.pfull*100,1), [('pentad', range(1,73)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
    wdudp_av = w_av*dudp_av

    duveddy_av = xr.DataArray( cfd(uved_av.values ,uved_av.lat*np.pi/180,2),   [('pentad', range(1,73)  ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat)])
    duveddy_av = duveddy_av/coslat/coslat/a
    duweddp_av = xr.DataArray( cfd(uwed_av.values,uwed_av.pfull*100,1),   [('pentad', range(1,73) ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat)])




    #evaluate sums of the terms
    mom_eddy = -(duveddy_av + duweddp_av)
    mom_mean = -(vdudy_av + wdudp_av)
    mom_sum = fv_av + ddamp_av + mom_eddy + mom_mean

    data_out = xr.Dataset({'fv_av': fv_av, 'ddamp_av': ddamp_av, 
                     'mom_eddy': mom_eddy, 'mom_mean': mom_mean, 'mom_sum': mom_sum , 'u_av': u_av, 'v_av': v_av   })


    return data_out



def mombudg_lon_pd_fn(run_fol,years, trange):
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
    v =     xr.DataArray(np.zeros((trl,2,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull[2:10:7]), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    w =     xr.DataArray(np.zeros((trl,17,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lon', rundata.lon), ('year', years )])
    uu =    xr.DataArray(np.zeros((trl,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ),  ('lon', rundata.lon), ('year', years )])
    uv =    xr.DataArray(np.zeros((trl,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uw =    xr.DataArray(np.zeros((trl,17,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull),  ('lon', rundata.lon), ('year', years )])
    phi =   xr.DataArray(np.zeros((trl,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ),  ('lon', rundata.lon), ('year', years )])
    ddamp = xr.DataArray(np.zeros((trl,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ),  ('lon', rundata.lon), ('year', years )])

    for year in years:
        print year
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) - 71
        u[:,:,:,:,year-years[0]]     = rundata.ucomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        v[:,:,:,:,year-years[0]]     = rundata.vcomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],2:10:7,:,:]
        w[:,:,:,year-years[0]]     = rundata.omega.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,37,:]
        uu[:,:,year-years[0]]    = rundata.ucomp_sq.groupby('pentad').mean(('time'))[trange[0]:trange[1],9,37,:]
        uv[:,:,:,year-years[0]]    = rundata.ucomp_vcomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],9,:,:]
        uw[:,:,:,year-years[0]]    = rundata.ucomp_omega.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,37,:]
        phi[:,:,year-years[0]]   = rundata.height.groupby('pentad').mean(('time'))[trange[0]:trange[1],9,37,:]
        ddamp[:,:,year-years[0]] = rundata.dt_ug_diffusion.groupby('pentad').mean(('time'))[trange[0]:trange[1],9,37,:]

    #Take averages: u, v, uv, w, phi, damping terms
    u_av = u.mean(('year'))
    v_av = v.mean(('year'))
    w_av = w.mean(('year'))
    uu_av = uu.mean(('year'))
    uv_av = uv.mean(('year'))
    uw_av = uw.mean(('year'))
    phi_av = phi.mean(('year'))*9.8

    ddamp_av = ddamp.mean(('year'))
    fv_av = v_av[:,1,37,:]*f[37]

    #Evaluate mean eddy products
    uued_av  =  uu_av - u_av[:,9,37,:]*u_av[:,9,37,:]
    uved_av  = (uv_av - u_av[:,9,:,:]*v_av[:,1,:,:])*coslat*coslat
    uwed_av  =  uw_av - u_av[:,:,37,:]*w_av
    

    #Evaluate gradients needed        
    dudx_av = xr.DataArray( cfd(u_av[:,9,37,:].values,u_av.lon*np.pi/180,1),   [('pentad', range(trange[0]+1,trange[1]+1) ), ('lon', u_av.lon )])
    ududx_av = u_av[:,9,37,:]*dudx_av/coslat[37]/a
    
    dudy_av = xr.DataArray( cfd((u_av[:,9,:,:]*coslat).values,u_av.lat*np.pi/180,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('lat', u_av.lat), ('lon', u_av.lon )])
    vdudy_av = v_av[:,1,37,:]*dudy_av[:,37,:]/coslat[37]/a
    
    dudp_av = xr.DataArray( cfd(u_av[:,:,37,:].values,u_av.pfull*100,1), [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lon', u_av.lon )])
    wdudp_av = w_av[:,9,:]*dudp_av[:,9,:]

    duueddx_av = xr.DataArray( cfd(uued_av.values,uued_av.lon*np.pi/180,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ),  ('lon', uued_av.lon )])
    duueddx_av = duueddx_av/coslat[37]/a
    duveddy_av = xr.DataArray( cfd(uved_av.values ,uved_av.lat*np.pi/180,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('lat', uved_av.lat), ('lon', uved_av.lon )])
    duveddy_av = duveddy_av[:,37,:]/coslat[37]/coslat[37]/a
    duweddp_av = xr.DataArray( cfd(uwed_av.values,uwed_av.pfull*100,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', uwed_av.pfull ), ('lon', uwed_av.lon )])
    duweddp_av = duweddp_av[:,9,:]
    
    dphidx_av = xr.DataArray( cfd(phi_av.values,phi_av.lon*np.pi/180,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('lon', phi_av.lon )])
    dphidx_av = -1*dphidx_av/coslat[37]/a


    #evaluate sums of the terms
    mom_eddy = -(duueddx_av + duveddy_av + duweddp_av)
    mom_mean = -(ududx_av + vdudy_av + wdudp_av)
    mom_sum = fv_av + ddamp_av + dphidx_av + mom_eddy + mom_mean

    data_out = xr.Dataset({'dphidx_av': dphidx_av, 'fv_av': fv_av, 'ddamp_av': ddamp_av, 
                     'mom_eddy': mom_eddy, 'mom_mean': mom_mean, 'mom_sum': mom_sum , 'u_av': u_av[:,9,37,:], 'v_av': v_av[:,1,37,:]   })
    wspd = xr.Dataset({'u_av': u_av[:,2,37,:], 'v_av': v_av[:,0,37,:]})

    return data_out, wspd





def mombudg_lev_pd_fn(run_fol,years, trange):
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
    v =     xr.DataArray(np.zeros((trl,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    w =     xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uu =    xr.DataArray(np.zeros((trl,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('lat', rundata.lat),  ('lon', rundata.lon), ('year', years )])
    uv =    xr.DataArray(np.zeros((trl,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    uw =    xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat),  ('lon', rundata.lon), ('year', years )])
    phi =   xr.DataArray(np.zeros((trl,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('lat', rundata.lat),  ('lon', rundata.lon), ('year', years )])
    ddamp = xr.DataArray(np.zeros((trl,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('lat', rundata.lat),  ('lon', rundata.lon), ('year', years )])

    for year in years:
        print year
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) - 71
        u[:,:,:,:,year-years[0]]     = rundata.ucomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        v[:,:,:,year-years[0]]     = rundata.vcomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],9,:,:]
        w[:,:,:,:,year-years[0]]     = rundata.omega.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        uu[:,:,:,year-years[0]]    = rundata.ucomp_sq.groupby('pentad').mean(('time'))[trange[0]:trange[1],9,:,:]
        uv[:,:,:,year-years[0]]    = rundata.ucomp_vcomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],9,:,:]
        uw[:,:,:,:,year-years[0]]    = rundata.ucomp_omega.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        phi[:,:,:,year-years[0]]   = rundata.height.groupby('pentad').mean(('time'))[trange[0]:trange[1],9,:,:]
        ddamp[:,:,:,year-years[0]] = rundata.dt_ug_diffusion.groupby('pentad').mean(('time'))[trange[0]:trange[1],9,:,:]

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
    uued_av  =  uu_av - u_av[:,9,:,:]*u_av[:,9,:,:]
    uved_av  = (uv_av - u_av[:,9,:,:]*v_av)*coslat*coslat
    uwed_av  =  uw_av - u_av*w_av
    

    #Evaluate gradients needed        
    dudx_av = xr.DataArray( cfd(u_av[:,9,:,:].values,u_av.lon*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1) ), ('lat', u_av.lat), ('lon', u_av.lon )])
    ududx_av = u_av[:,9,:,:]*dudx_av/coslat/a
    
    dudy_av = xr.DataArray( cfd((u_av[:,9,:,:]*coslat).values,u_av.lat*np.pi/180,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('lat', u_av.lat), ('lon', u_av.lon )])
    vdudy_av = v_av*dudy_av/coslat/a
    
    dudp_av = xr.DataArray( cfd(u_av.values,u_av.pfull*100,1), [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    wdudp_av = w_av[:,9,:,:]*dudp_av[:,9,:,:]

    duueddx_av = xr.DataArray( cfd(uued_av.values,uued_av.lon*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('lat', u_av.lat),  ('lon', uued_av.lon )])
    duueddx_av = duueddx_av/coslat/a
    duveddy_av = xr.DataArray( cfd(uved_av.values ,uved_av.lat*np.pi/180,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('lat', uved_av.lat), ('lon', uved_av.lon )])
    duveddy_av = duveddy_av/coslat/coslat/a
    duweddp_av = xr.DataArray( cfd(uwed_av.values,uwed_av.pfull*100,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', uwed_av.pfull ), ('lat', u_av.lat), ('lon', uwed_av.lon )])
    duweddp_av = duweddp_av[:,9,:,:]
    
    dphidx_av = xr.DataArray( cfd(phi_av.values,phi_av.lon*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('lat', u_av.lat), ('lon', phi_av.lon )])
    dphidx_av = -1*dphidx_av/coslat/a


    #evaluate sums of the terms
    mom_eddy = -(duueddx_av + duveddy_av + duweddp_av)
    mom_mean = -(ududx_av + vdudy_av + wdudp_av)
    mom_sum = fv_av + ddamp_av + dphidx_av + mom_eddy + mom_mean
    
    data_out = xr.Dataset({'dphidx_av': dphidx_av, 'fv_av': fv_av, 'ddamp_av': ddamp_av, 
                     'mom_eddy': mom_eddy, 'mom_mean': mom_mean, 'mom_sum': mom_sum , 'u_av': u_av[:,9,:,:], 'v_av': v_av   })


    return data_out







   
if __name__ == "__main__":

        #set run name
        inp_fol = 'flat_10m'
        years = [21]
	mombudg_lev_pd_fn(inp_fol,years,[0,72])
