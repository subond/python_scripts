#Load in u, v, omega, geopotential, evaluate the terms in the zonal momentum budget
#NB may also need to account for damping terms

#Start by doing this for an annual average, can extend after

import sys
sys.path.append('/scratch/rg419/workdir_moist/python_bin')
from time_av_xr import load_year_xr
import numpy as np
import matplotlib.pyplot as plt
from finite_difference import cfd
import xarray as xr

#set run name

def vmombudg_fn(inp_fol,nproc,dim):

#inp_fol = experiment name
#nproc   = no of processors
#dim     = dimension to take eddies around (time or lon)

    #load in run data
    run_fol = inp_fol + '/np' + nproc
    year = 5
    rundata = load_year_xr(run_fol, year, pinterp=True)

    #Define constants
    omega = 7.2921150e-5
    a= 6376.0e3 #radius used in model
    coslat = np.cos(rundata.lat * np.pi/180)
    f = 2 * omega * np.sin(rundata.lat *np.pi/180)

    #Take first averages: u, v, uv, w, phi, damping terms
    u_av = rundata.ucomp.mean((dim))
    v_av = rundata.vcomp.mean((dim))
    w_av = rundata.omega.mean((dim))
    phi_av = rundata.height.mean((dim))*9.8

    #terms that don't need differentiating etc can be immediately averaged
    ddamp_tzav = rundata.dt_vg_diffusion.mean(('time','lon'))
    rdamp_tzav = rundata.vdt_rdamp.mean(('time','lon'))
    fu_tzav = f*rundata.ucomp.mean(('time','lon'))

    #look at eddies from chosen average
    u_ed  = rundata.ucomp - u_av
    v_ed  = rundata.vcomp - v_av
    w_ed  = rundata.omega - w_av

    #calculate eddy products
    vved_av = (v_ed*v_ed).mean((dim))*coslat*coslat
    uved_av = (u_ed*v_ed).mean((dim))
    vwed_av = (v_ed*w_ed).mean((dim))

    #Evaluate gradients needed
    if dim == 'time':
        dvdx_av = xr.DataArray( cfd(v_av.values,u_av.lon*np.pi/180,2),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        udvdx_av = u_av*dvdx_av/coslat/a
        dvdy_av = xr.DataArray( cfd((v_av*coslat).values,u_av.lat*np.pi/180,1),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        vdvdy_av = v_av*dvdy_av/coslat/a
        dvdp_av = xr.DataArray( cfd(v_av.values,u_av.pfull*100,0), [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        wdvdp_av = w_av*dvdp_av

        duveddx_av = xr.DataArray( cfd(uved_av.values,u_av.lon*np.pi/180,2),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        duveddx_av = duveddx_av/coslat/a
        dvveddy_av = xr.DataArray( cfd(vved_av.values ,u_av.lat*np.pi/180,1),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
        dvveddy_av = dvveddy_av/coslat/coslat/a
        dvweddp_av = xr.DataArray( cfd(vwed_av.values,u_av.pfull*100,0),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])

        dphidy_av = xr.DataArray( cfd(phi_av.values,phi_av.lat*np.pi/180,1),   [('pfull', phi_av.pfull ), ('lat', phi_av.lat), ('lon', phi_av.lon )])
        dphidy_av = dphidy_av/a

        dphidy_tzav = dphidy_av.mean(('lon'))

        udvdx_tzav = udvdx_av.mean(('lon'))
        vdvdy_tzav = vdvdy_av.mean(('lon'))
        wdvdp_tzav = wdvdp_av.mean(('lon'))

        duveddx_tzav = duveddx_av.mean(('lon'))
        dvveddy_tzav = dvveddy_av.mean(('lon'))
        dvweddp_tzav = dvweddp_av.mean(('lon'))

    #evaluate a sum of the terms
    mom_sum = dphidy_tzav + fu_tzav + udvdx_tzav + vdvdy_tzav + wdvdp_tzav + duveddx_tzav + dvveddy_tzav + dvweddp_tzav - ddamp_tzav
    

    #Plot up terms to check
    plt.figure(1)
    (dphidy_tzav+fu_tzav).plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(2)
    udvdx_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(3)
    vdvdy_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(4)
    wdvdp_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(5)
    duveddx_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(6)
    dvveddy_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(7)
    dvweddp_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(8)
    mom_sum.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(9)
    ddamp_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)

    #produce Dima&Wallace plots
    #plt.figure(1)
    #(-duveddy_tzav).plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    #plt.figure(2)
    #(fv_tzav-vdudy_tzav).plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    #plt.figure(3)
    #(fv_tzav-vdudy_tzav-wdudp_tzav-duveddy_tzav-duweddp_tzav).plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.show()

    return 

   
if __name__ == "__main__":

        #set run name
        inp_fol = 'aquaplanet_10m_diags'
        nproc   = '16'

	vmombudg_fn(inp_fol,nproc,'time')
