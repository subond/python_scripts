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

def mombudg_fn(inp_fol,nproc,dim):

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
    uu_av = rundata.ucomp_sq.mean((dim))
    uv_av = rundata.ucomp_vcomp.mean((dim))
    uw_av = rundata.ucomp_omega.mean((dim))
    phi_av = rundata.height.mean((dim))*9.8

    #terms that don't need differentiating etc can be immediately averaged
    ddamp_tzav = rundata.dt_ug_diffusion.mean(('time','lon'))
    rdamp_tzav = rundata.udt_rdamp.mean(('time','lon'))
    fv_tzav = f*rundata.vcomp.mean(('time','lon'))

    #calculate eddy products
    uued_av  =  uu_av - u_av*u_av
    uved_av  = (uv_av - u_av*v_av)*coslat*coslat
    uwed_av  =  uw_av - u_av*w_av

    #Evaluate gradients needed
    if dim == 'time':
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
        dphidx_av = dphidx_av/coslat/a

        dphidx_tzav = dphidx_av.mean(('lon'))

        ududx_tzav = ududx_av.mean(('lon'))
        vdudy_tzav = vdudy_av.mean(('lon'))
        wdudp_tzav = wdudp_av.mean(('lon'))

        duueddx_tzav = duueddx_av.mean(('lon'))
        duveddy_tzav = duveddy_av.mean(('lon'))
        duweddp_tzav = duweddp_av.mean(('lon'))

    elif dim == 'lon':
        dudy_av = xr.DataArray( cfd((u_av*coslat).values,u_av.lat*np.pi/180,2),   [('time', u_av.time ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
        vdudy_av = v_av*dudy_av/coslat/a
        dudp_av = xr.DataArray( cfd(u_av.values,u_av.pfull*100,1), [('time', u_av.time ), ('pfull', u_av.pfull ), ('lat', u_av.lat)])
        wdudp_av = w_av*dudp_av

        duveddy_av = xr.DataArray( cfd(uved_av.values ,uved_av.lat*np.pi/180,2),   [('time', u_av.time ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat)])
        duveddy_av = duveddy_av/coslat/coslat/a
        duweddp_av = xr.DataArray( cfd(uwed_av.values,uwed_av.pfull*100,1),   [('time', u_av.time ), ('pfull', uued_av.pfull ), ('lat', uued_av.lat)])

        vdudy_tzav = vdudy_av.mean(('time'))
        wdudp_tzav = wdudp_av.mean(('time'))
        duveddy_tzav = duveddy_av.mean(('time'))
        duweddp_tzav = duweddp_av.mean(('time'))

        dphidx_tzav =  xr.DataArray( np.zeros((40,64)),   [('pfull', u_av.pfull ), ('lat', u_av.lat)])
        ududx_tzav =   xr.DataArray( np.zeros((40,64)),   [('pfull', u_av.pfull ), ('lat', u_av.lat)])
        duueddx_tzav = xr.DataArray( np.zeros((40,64)),   [('pfull', u_av.pfull ), ('lat', u_av.lat)])

    #evaluate a sum of the terms
    mom_sum = dphidx_tzav - fv_tzav + ududx_tzav + vdudy_tzav + wdudp_tzav + duueddx_tzav + duveddy_tzav + duweddp_tzav - ddamp_tzav

    #produce output dataset
    data_out = xr.Dataset({'dphidx_tzav': dphidx_tzav, 'fv_tzav': fv_tzav, 'ddamp_tzav': ddamp_tzav, 'rdamp_tzav': rdamp_tzav,
                     'ududx_tzav': ududx_tzav, 'vdudy_tzav': vdudy_tzav, 'wdudp_tzav': wdudp_tzav,
                     'duueddx_tzav': duueddx_tzav, 'duveddy_tzav': duveddy_tzav, 'duweddp_tzav': duweddp_tzav  })
    

    #Plot up terms to check
    plt.figure(1)
    dphidx_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(2)
    fv_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(3)
    ududx_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(4)
    vdudy_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(5)
    wdudp_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(6)
    duueddx_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(7)
    duveddy_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(8)
    duweddp_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(9)
    mom_sum.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.figure(10)
    ddamp_tzav.plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)

    #produce Dima&Wallace plots
    #plt.figure(1)
    #(-duveddy_tzav).plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    #plt.figure(2)
    #(fv_tzav-vdudy_tzav).plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    #plt.figure(3)
    #(fv_tzav-vdudy_tzav-wdudp_tzav-duveddy_tzav-duweddp_tzav).plot.contourf(x='lat', y='pfull',levels=np.arange(-0.0001,0.0001,0.00001), yincrease=False)
    plt.show()

    return data_out


#function to evaluate vertically integrated terms
def mombudg_vint_fn(inp_fol,nproc,dim):
    rundata, phalf = mombudg_fn(inp_fol,nproc,dim)
    pdiffs = xr.DataArray(phalf.diff('phalf').values, [('pfull', rundata.pfull )])

    fv_vint = (rundata.fv_tzav * pdiffs*100).sum(('pfull'))/9.8
    ddamp_vint = (rundata.ddamp_tzav * pdiffs*100).sum(('pfull'))/9.8
    rdamp_vint = (rundata.rdamp_tzav * pdiffs*100).sum(('pfull'))/9.8
    ududx_vint = (rundata.ududx_tzav * pdiffs*100).sum(('pfull'))/9.8
    vdudy_vint = (rundata.vdudy_tzav * pdiffs*100).sum(('pfull'))/9.8
    wdudp_vint = (rundata.wdudp_tzav * pdiffs*100).sum(('pfull'))/9.8
    duueddx_vint = (rundata.duueddx_tzav * pdiffs*100).sum(('pfull'))/9.8
    duveddy_vint = (rundata.duveddy_tzav * pdiffs*100).sum(('pfull'))/9.8
    duweddp_vint = (rundata.duweddp_tzav * pdiffs*100).sum(('pfull'))/9.8

    mom_sum = ududx_vint + vdudy_vint + wdudp_vint + duueddx_vint + duveddy_vint + duweddp_vint - ddamp_vint - rdamp_vint

    fv_vint.plot()
    #rdamp_vint.plot()
    ddamp_vint.plot()
    (vdudy_vint+wdudp_vint).plot()
    #wdudp_vint.plot()
    duveddy_vint.plot()
    duweddp_vint.plot()
    mom_sum.plot()
    plt.legend(['fv','ddamp','C_zonal','duvdy','duwdp','sum'])
    plt.show()
   
if __name__ == "__main__":

        #set run name
        inp_fol = 'aquaplanet_10m_diags'
        nproc   = '16'

	mombudg_fn(inp_fol,nproc,'time')
