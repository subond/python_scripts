# Will be a re-write of thetaint_version4_2015_earth.pro.
import numpy as np
import interpolation as interp
import analyse as ana


## ------------------------------------- These variables are inputs---------------------------->|These variables are outputs-->|
#PRO thetaint_version4_2015_earth, u,v,ps,temp,thetapre,lonar,latar,sigar,lonartot,latartot,PV,utheta,vtheta,vort,dvdx,dudy

def pv_calc(thd_data,vars_list,thetapre,lonar,latar,pfull,lonartot,latartot, planet_params):

#planet_params=[rm,g,rct,rctg,pref,kappa,omega]

	rm=planet_params[0]
	g=planet_params[1]
	rct=planet_params[2]
	rctg=planet_params[3]
	pref=planet_params[4]
	kappa=planet_params[5]
	omega=planet_params[6]

	nlons=len(lonar)
	nlats=len(latar)
	nlevs=len(pfull)
	ntheta=len(thetapre)
	ntime_tot=len(np.squeeze(thd_data[0,:,0,0,0]))
	
	n_arr=[nlons,nlats,nlevs,ntheta,ntime_tot]

	temp_index=vars_list.index('temp')
	u_index=vars_list.index('ucomp')
	v_index=vars_list.index('vcomp')
	ptemp_index=vars_list.index('ptemp')

	temp=np.squeeze(thd_data[temp_index,:,:,:,:])
	ucomp=np.squeeze(thd_data[u_index,:,:,:,:])
	vcomp=np.squeeze(thd_data[v_index,:,:,:,:])
	ptemp=np.squeeze(thd_data[ptemp_index,:,:,:,:])
	
	thd_temp=np.zeros((3,ntime_tot,nlevs,nlats,nlons))


#*****************************************************************************************************************************************************************************************************
# Description of criteria for longitude arrays and latitude arrays IF they are a SUBSET of the model's global coverage. (i.e. if the model goes from -87.5, 87.5 latitude, but you only pass data between -30,30, etc.)
# If subroutine gets passed a longitude array that is a subset of the true global dataset in longitude, the PV will not be calculated at the longitudes at the edge of the array, as it cannot calculate the relative vorticity using the central differencing method used in what follows. I.e. if you pass the routine an array of data in longitude corresponding to:
#
# lonar=[-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.]
#
# then PV will only be caluclated at longitudes:
# [-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.]
#
# So, if you wish to have PV calculated at +/-25 in this example, you'll need to provide data in the longitude range:
# lonar=[-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.]
#
#If, however, you pass data that is periodic in longitude, the PV will be calculated at all longitudes.
# 
# As the latitudinal direction is not periodic in the same way the longitudinal direction is, the PV will ALWAYS be outputted on an array that is smaller by 2 points in the latitudinal direction, even if the global dataset is passed. I.e. if
# latar=[-87.5,-85.,-82.5,......82.5,85.,87.5,]
# then the PV will only be outputted at:
# [-85.,-82.5,......82.5,85.]
#
# If a subset is passed in latitude, the same procedure is applied as is described for longitude above
# *****************************************************************************************************************************************************************************************************



	#Need to see if prescribed theta levels are equally spaced.

	dtheta_test=thetapre[1:4]-thetapre[0:3]
	theta_test_result=all(np.absolute(x-dtheta_test[0]) < 1e-10 for x in dtheta_test)

	if not theta_test_result:
		print 'Cannot currently cope with non-equally-spaced theta levels as input. ERROR.'
		return

#--------------------------------This is test to see if I am dealing with the entire global dataset, or a subset of it--------------------------------
	correctlatstart=0	
	correctlatend=0
	wholelat=1

	if len(latar) != len(latartot):
		wholelat=0
		if latar[0] != np.max(latartot):
			correctlatstart=1
		if latar[-1] != np.min(latartot):
			correctlatstart=1

	correctlonstart=0	
	correctlonend=0
	wholelon=1

	if len(lonar) != len(lonartot):
		wholelon=0
		if lonar[0] != np.max(lonartot):
			correctlonstart=1
		if lonar[-1] != np.min(lonartot):
			correctlonstart=1
#--------------------------------------------------------------------------------------------------------------------------------------------------------

	#First thing is to interpolate the winds onto theta levels. We already have ptemp...

	print 'interpolating winds'

	utheta=np.zeros((ntime_tot,ntheta,nlats,nlons))
	vtheta=np.zeros((ntime_tot,ntheta,nlats,nlons))
	temptheta=np.zeros((ntime_tot,ntheta,nlats,nlons))


	in_u=np.zeros(ntheta)
	in_v=np.zeros(ntheta)
	in_t=np.zeros(ntheta)
	thetalev=np.zeros(ntheta)

	time=range(ntime_tot)
	lats=range(nlats)
	lons=range(nlons)
	levs=range(nlevs)
	thetalevs=range(ntheta)


	for nt in time:
		for nlat in lats:
			for nlon in lons:

				tprof=temp[nt,:,nlat,nlon]
				uprof=ucomp[nt,:,nlat,nlon]
				vprof=vcomp[nt,:,nlat,nlon]
				thetalev=ptemp[nt,:,nlat,nlon]


				for nlev in thetalevs:
					lolev, vw = interp.wei_vert(0,thetalev,thetapre[nlev])

					if lolev >= 0:
						in_t[nlev]=tprof[lolev]+(tprof[lolev+1]-tprof[lolev])*vw
						in_u[nlev]=uprof[lolev]+(uprof[lolev+1]-uprof[lolev])*vw
						in_v[nlev]=vprof[lolev]+(vprof[lolev+1]-vprof[lolev])*vw

					else:
						in_t[nlev]=0
						in_u[nlev]=0
						in_v[nlev]=0

					temptheta[nt,:,nlat,nlon]=in_t
					utheta[nt,:,nlat,nlon]=in_u
					vtheta[nt,:,nlat,nlon]=in_v

		print nt


# 	vars_list.append('utheta')
# 	vars_list.append('vtheta')
# 	vars_list.append('temptheta')

# 	thd_temp=np.append(thd_data,[utheta],axis=0)
# 	thd_temp=np.append(thd_temp,[vtheta],axis=0)
# 	thd_temp=np.append(thd_temp,[temptheta],axis=0)

	vars_list_theta=['utheta','vtheta','temptheta']
	 
	thd_temp[0,:,:,:,:]=utheta
	thd_temp[1,:,:,:,:]=vtheta
	thd_temp[2,:,:,:,:]=temptheta
	
	
	
# --------------------------------------Now have winds on Theta Levels-----------------------------------------------------
# ---------------------------- Now want to calculate Vorticity on Theta Levels---------------------------------------------


	vort_theta=ana.rel_vort_calc(utheta,vtheta,lonar,latar,rm,n_arr,wholelat,wholelon)

	vars_list_theta.append('vort_theta')

	thd_temp=np.append(thd_temp,[vort_theta],axis=0)

# Why not just subset thd_temp and others whilst you're here? Might make it easier than doing the below!



# ; ******************************************************Begin enforced sub-setting of y data**************************************************************
# ; At this point, because we have not calculated the y gradient of u at the two edge latitudinal gridpoints, and therefore we cannot claim to have a vorticity, or a PV defined there. So, we throw all of the data at the two final gridpoints in latitude away. Do that with the whole dataset, so as to make the rest of the routine think it is simply dealing with a dataset of nlat-2 gridpoints in y, as opposed to one with nlat.
# ; stop
# 
# 
	print thd_temp.shape
	print thd_data.shape
	print latar.shape

	thd_temp=thd_temp[:,:,:,1:-2,:]
	thd_data=thd_data[:,:,:,1:-2,:]
	latar=latar[1:-2]

	print thd_temp.shape
	print thd_data.shape
	print latar.shape
	
	#s This array is too small - why?

# ; ******************************************************End enforced sub-setting of y data**************************************************************
# 
# 
# ; ******************************************************Begin optional sub-setting of x data**************************************************************
# ; If we are given the whole longitudinal extent of the data (i.e. data that covers the whole globe in the longitudinal direction), then x gradients at the two 'edges' of the domain (i.e. lonar(0) and lonar(-1)), can be calculated using this periodicity in x. So, there is no enforced sub-setting of the data in the x direction, as there is in y. However, if the whole longitudinal dataset is not passed to the routine (i.e. wholelon=0) then the x gradient of the data cannot be calculated at the extrema in x, and so the data must be sub-sampled in x.
# 


	if wholelon !=1:
		thd_temp=thd_temp[:,:,:,:,1:-2]
		thd_data=thd_data[:,:,:,:,1:-2]
		lonar=lonar[1:-2]

# ; ******************************************************End optional sub-setting of x data**************************************************************
# 



			
	
	return vars_list_theta, thd_temp, thd_data, vars_list, lonar, latar
