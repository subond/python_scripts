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

def mombudg_2d_fn(inp_fol,nproc):

#inp_fol = experiment name
#nproc   = no of processors

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
    u_av = rundata.ucomp.mean(('time'))
    v_av = rundata.vcomp.mean(('time'))
    w_av = rundata.omega.mean(('time'))
    phi_av = rundata.height.mean(('time'))*9.8

    #terms that don't need differentiating etc can be immediately averaged
    ddamp_av = rundata.dt_vg_diffusion.mean(('time'))
    fu_av = f*rundata.ucomp.mean(('time'))

    #look at eddies from chosen average
    u_ed  = rundata.ucomp - u_av
    v_ed  = rundata.vcomp - v_av
    w_ed  = rundata.omega - w_av

    #calculate eddy products
    uved_av = (u_ed*v_ed).mean(('time'))
    vved_av = (v_ed*v_ed).mean(('time'))*coslat*coslat
    vwed_av = (v_ed*w_ed).mean(('time'))

    #Evaluate gradients needed
    dvdx_av = xr.DataArray( cfd(v_av.values,u_av.lon*np.pi/180,2),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    udvdx_av = u_av*dvdx_av/coslat/a
    dvdy_av = xr.DataArray( cfd((v_av*coslat).values,u_av.lat*np.pi/180,1),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    vdvdy_av = v_av*dvdy_av/coslat/a
    dvdp_av = xr.DataArray( cfd(u_av.values,u_av.pfull*100,0), [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    wdvdp_av = w_av*dvdp_av

    duveddx_av = xr.DataArray( cfd(uved_av.values,u_av.lon*np.pi/180,2),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    duveddx_av = duveddx_av/coslat/a
    dvveddy_av = xr.DataArray( cfd(vved_av.values ,u_av.lat*np.pi/180,1),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])
    dvveddy_av = dvveddy_av/coslat/coslat/a
    dvweddp_av = xr.DataArray( cfd(vwed_av.values,u_av.pfull*100,0),   [('pfull', u_av.pfull ), ('lat', u_av.lat), ('lon', u_av.lon )])

    dphidy_av = xr.DataArray( cfd(phi_av.values,phi_av.lat*np.pi/180,1),   [('pfull', phi_av.pfull ), ('lat', phi_av.lat), ('lon', phi_av.lon )])
    dphidy_av = dphidy_av/a


    #evaluate sums of the terms
    mom_eddy = -(duveddx_av + dvveddy_av + dvweddp_av)
    mom_mean = -(udvdx_av + vdvdy_av + wdvdp_av)
    mom_sum = - fu_av + ddamp_av - dphidy_av + mom_eddy + mom_mean

    #Make a test plot to see if local balance appears to be achieved...
    plt.figure(1)
    (dphidy_av[9,:,:]+fu_av[:,9,:]).plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    plt.ylim(-30,30)
    plt.figure(2)
    fu_av[:,9,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    plt.ylim(-30,30)
    plt.figure(3)
    mom_eddy[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    plt.ylim(-30,30)
    plt.figure(4)
    mom_mean[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    plt.ylim(-30,30)
    plt.figure(5)
    ddamp_av[9,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    plt.ylim(-30,30)
    plt.figure(6)
    mom_sum[:,9,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0002,0.0002,0.00001))
    plt.ylim(-30,30)

    #plt.figure(7)
    #(dphidy_av[2,:,:]+fu_av[:,2,:]).plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.0001,0.00001))
    #plt.ylim(-30,30)
    #plt.figure(8)
    #fu_av[:,2,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.0001,0.00001))
    #plt.ylim(-30,30)
    #plt.figure(9)
    #mom_eddy[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.0001,0.00001))
    #plt.ylim(-30,30)
    #plt.figure(10)
    #mom_mean[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.0001,0.00001))
    #plt.ylim(-30,30)
    #plt.figure(11)
    #ddamp_av[2,:,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.0001,0.00001))
    #plt.ylim(-30,30)
    #plt.figure(12)
    #mom_sum[:,2,:].plot.contourf(x='lon', y='lat',levels=np.arange(-0.0001,0.0001,0.00001))
    #plt.ylim(-30,30)
    plt.show()

    return 


   
if __name__ == "__main__":

        #set run name
        inp_fol = 'aquaplanet_10m_diags'
        nproc   = '16'

	mombudg_2d_fn(inp_fol,nproc)
