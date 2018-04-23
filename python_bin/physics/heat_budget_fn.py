#Load in u, v, omega, geopotential, evaluate the terms in the zonal momentum budget

from data_handling import load_year_xr
import numpy as np
import matplotlib.pyplot as plt
from finite_difference import cfd
import xarray as xr

def heatbudg_2d_pd_fn(run_fol,years, trange, heatbudg=True, humbudg=False):
    year = years[0]
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['pentad'] = (rundata.time // 5) - 71
    mngrp = rundata.ucomp.groupby('pentad').mean(('time'))
    trl = trange[1]-trange[0]
    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(rundata.lat * np.pi/180)

    #Initialise arrays to load into
    u =     xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    v =     xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    w =     xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    if heatbudg == True:
        t =    xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    if humbudg == True:
        q =    xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    if heatbudg==False and humbudg==False:
        print 'No budget specified'
        return
        
    for year in years:
        print year
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) - 71
        u[:,:,:,:,year-years[0]]     = rundata.ucomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        v[:,:,:,:,year-years[0]]     = rundata.vcomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        w[:,:,:,:,year-years[0]]     = rundata.omega.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        if heatbudg==True:
            t[:,:,:,:,year-years[0]]     = rundata.temp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        if humbudg==True:
            q[:,:,:,:,year-years[0]]     = rundata.sphum.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]

    #Take averages: u, v, uv, w, phi, damping terms
    u_av = u.mean(('year'))
    v_av = v.mean(('year'))
    w_av = w.mean(('year'))

    #Evaluate eddy quantities
    u_ed = u - u_av
    v_ed = v - v_av
    w_ed = w - w_av
    
    if heatbudg == True:
        t_av = t.mean(('year'))
        t_ed = t - t_av
        ut_ed  =  (u_ed*t_ed).mean(('year'))
        vt_ed  =  (v_ed*t_ed).mean(('year'))*coslat
        wt_ed  =  (w_ed*t_ed).mean(('year'))
        
        #Evaluate gradients needed
        dtdx_av = xr.DataArray( cfd(t_av.values,t_av.lon*np.pi/180,3),   [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        udtdx_av = u_av*dtdx_av/coslat/a
        dtdy_av = xr.DataArray( cfd(t_av.values,t_av.lat*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        vdtdy_av = v_av*dtdy_av/a
        dtdp_av = xr.DataArray( cfd(t_av.values,t_av.pfull*100,1), [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        wdtdp_av = w_av*dtdp_av

        duteddx_av = xr.DataArray( cfd(ut_ed.values,ut_ed.lon*np.pi/180,3),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        duteddx_av = duteddx_av/coslat/a
        dvteddy_av = xr.DataArray( cfd(vt_ed.values ,vt_ed.lat*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        dvteddy_av = dvteddy_av/coslat/a
        dwteddp_av = xr.DataArray( cfd(wt_ed.values,wt_ed.pfull*100,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        
        #evaluate sums of the terms
        heat_eddy = -(duteddx_av + dvteddy_av + dwteddp_av)
        heat_mean = -(udtdx_av + vdtdy_av + wdtdp_av)
        heat_sum = heat_eddy + heat_mean
        
    if humbudg == True:
        q_av = q.mean(('year'))
        q_ed = q - q_av
        uq_ed  =  (u_ed*q_ed).mean(('year'))
        vq_ed  =  (v_ed*q_ed).mean(('year'))*coslat
        wq_ed  =  (w_ed*q_ed).mean(('year'))
        
        #Evaluate gradients needed
        dqdx_av = xr.DataArray( cfd(q_av.values,q_av.lon*np.pi/180,3),   [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        udqdx_av = u_av*dqdx_av/coslat/a
        dqdy_av = xr.DataArray( cfd(q_av.values,q_av.lat*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        vdqdy_av = v_av*dqdy_av/a
        dqdp_av = xr.DataArray( cfd(q_av.values,q_av.pfull*100,1), [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        wdqdp_av = w_av*dqdp_av

        duqeddx_av = xr.DataArray( cfd(uq_ed.values,uq_ed.lon*np.pi/180,3),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        duqeddx_av = duqeddx_av/coslat/a
        dvqeddy_av = xr.DataArray( cfd(vq_ed.values,vq_ed.lat*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        dvqeddy_av = dvqeddy_av/coslat/a
        dwqeddp_av = xr.DataArray( cfd(wq_ed.values,wq_ed.pfull*100,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        
        #evaluate sums of the terms
        hum_eddy = -(duqeddx_av + dvqeddy_av + dwqeddp_av)
        hum_mean = -(udqdx_av + vdqdy_av + wdqdp_av)
        hum_sum = hum_eddy + hum_mean

    if heatbudg == True and humbudg == True:
        data_out = xr.Dataset({'hum_eddy': hum_eddy, 'hum_mean': hum_mean, 'hum_sum': hum_sum ,
                            'heat_eddy': heat_eddy, 'heat_mean': heat_mean, 'heat_sum': heat_sum , 'u_av': u_av, 'v_av': v_av   })
    elif heatbudg == True:
        data_out = xr.Dataset({'heat_eddy': heat_eddy, 'heat_mean': heat_mean, 'heat_sum': heat_sum , 'u_av': u_av, 'v_av': v_av   })
        
    elif humbudg == True:
        data_out = xr.Dataset({'hum_eddy': hum_eddy, 'hum_mean': hum_mean, 'hum_sum': hum_sum , 'u_av': u_av, 'v_av': v_av   })

    return data_out





def heatbudg_zon_pd_fn(run_fol,years, trange, heatbudg=True, humbudg=False):
    year = years[0]
    rundata = load_year_xr(run_fol, year, pinterp=True)
    rundata.coords['pentad'] = (rundata.time // 5) - 71
    mngrp = rundata.ucomp.groupby('pentad').mean(('time'))
    trl = trange[1]-trange[0]
    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(rundata.lat * np.pi/180)

    #Initialise arrays to load into
    u =     xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    v =     xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    w =     xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    if heatbudg == True:
        t =    xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    if humbudg == True:
        q =    xr.DataArray(np.zeros((trl,17,64,128,len(years))), [('pentad', range(trange[0]+1,trange[1]+1) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years )])
    if heatbudg==False and humbudg==False:
        print 'No budget specified'
        return
        
    for year in years:
        print year
        rundata = load_year_xr(run_fol, year, pinterp=True)
        rundata.coords['pentad'] = (rundata.time // 5) - 71
        u[:,:,:,:,year-years[0]]     = rundata.ucomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        v[:,:,:,:,year-years[0]]     = rundata.vcomp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        w[:,:,:,:,year-years[0]]     = rundata.omega.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        if heatbudg==True:
            t[:,:,:,:,year-years[0]]     = rundata.temp.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]
        if humbudg==True:
            q[:,:,:,:,year-years[0]]     = rundata.sphum.groupby('pentad').mean(('time'))[trange[0]:trange[1],:,:,:]

    #Take averages: u, v, uv, w, phi, damping terms
    u_av = u.mean(('year','lon'))
    v_av = v.mean(('year','lon'))
    w_av = w.mean(('year','lon'))

    #Evaluate eddy quantities
    u_ed = u - u_av
    v_ed = v - v_av
    w_ed = w - w_av
    
    if heatbudg == True:
        t_av = t.mean(('year','lon'))
        t_ed = t - t_av
        ut_ed  =  (u_ed*t_ed).mean(('year','lon'))
        vt_ed  =  (v_ed*t_ed).mean(('year','lon'))*coslat
        wt_ed  =  (w_ed*t_ed).mean(('year','lon'))
        
        #Evaluate gradients needed
        dtdy_av = xr.DataArray( cfd(t_av.values,t_av.lat*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
        vdtdy_av = v_av*dtdy_av/a
        dtdp_av = xr.DataArray( cfd(t_av.values,t_av.pfull*100,1), [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
        wdtdp_av = w_av*dtdp_av

        dvteddy_av = xr.DataArray( cfd(vt_ed.values ,vt_ed.lat*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
        dvteddy_av = dvteddy_av/coslat/a
        dwteddp_av = xr.DataArray( cfd(wt_ed.values,wt_ed.pfull*100,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
        
        #evaluate sums of the terms
        heat_eddy = -(dvteddy_av + dwteddp_av)
        heat_mean = -(vdtdy_av + wdtdp_av)
        heat_sum = heat_eddy + heat_mean
        
    if humbudg == True:
        q_av = q.mean(('year','lon'))
        q_ed = q - q_av
        uq_ed  =  (u_ed*q_ed).mean(('year','lon'))
        vq_ed  =  (v_ed*q_ed).mean(('year','lon'))*coslat
        wq_ed  =  (w_ed*q_ed).mean(('year','lon'))
        
        #Evaluate gradients needed
        dqdy_av = xr.DataArray( cfd(q_av.values,q_av.lat*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
        vdqdy_av = v_av*dqdy_av/a
        dqdp_av = xr.DataArray( cfd(q_av.values,q_av.pfull*100,1), [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
        wdqdp_av = w_av*dqdp_av

        dvqeddy_av = xr.DataArray( cfd(vq_ed.values,vq_ed.lat*np.pi/180,2),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
        dvqeddy_av = dvqeddy_av/coslat/a
        dwqeddp_av = xr.DataArray( cfd(wq_ed.values,wq_ed.pfull*100,1),   [('pentad', range(trange[0]+1,trange[1]+1)  ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
        
        #evaluate sums of the terms
        hum_eddy = -(dvqeddy_av + dwqeddp_av)
        hum_mean = -(vdqdy_av + wdqdp_av)
        hum_sum = hum_eddy + hum_mean

    if heatbudg == True and humbudg == True:
        data_out = xr.Dataset({'hum_eddy': hum_eddy, 'hum_mean': hum_mean, 'hum_sum': hum_sum ,
                            'heat_eddy': heat_eddy, 'heat_mean': heat_mean, 'heat_sum': heat_sum , 'u_av': u_av, 'v_av': v_av   })
    elif heatbudg == True:
        data_out = xr.Dataset({'heat_eddy': heat_eddy, 'heat_mean': heat_mean, 'heat_sum': heat_sum , 'u_av': u_av, 'v_av': v_av   })
        
    elif humbudg == True:
        data_out = xr.Dataset({'hum_eddy': hum_eddy, 'hum_mean': hum_mean, 'hum_sum': hum_sum , 'u_av': u_av, 'v_av': v_av   })

    return data_out





   
if __name__ == "__main__":

        #set run name
        inp_fol = 'aquaplanet_2m'
        years = range(21,22)
	heatbudg_2d_pd_fn(inp_fol,years,[19,61],humbudg=True)
