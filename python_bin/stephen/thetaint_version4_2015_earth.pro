; IDL PROCEDURE THETAINT
;
;-----DESCRIPTION
; This script calculates potential vorticity from 3D gridded data of winds and temperature provided on model sigma levels, and surface pressure.
; Potential temperatures (theta) are first calculated on the model sigma levels, then zonal and meridional winds are interpolatedd onto predefined potential temperature levels  
; Finite-difference spatial derivatives of winds are calculated, from which relative and potential vorticity fields on potential temperature surfaces are finally derived.
; The planet's geometry is assumed to be spherical.
; 
;-----COPYLEFT
; This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License (http://www.gnu.org/copyleft/gpl.html)
;
;-----HISTORY AND ACKNOWLEDGEMENTS
; IDL script originally written by Stephen I. Thomson in February 2011, with vertical interpolation adapted from an IDL script originally written by Luca Montabone.  
; Modifications and improvements by Stephen I. Thomson and Luca Montabone in July 2011.
; Further modifications and improvements by Stephen I. Thomson in June 2014.
; Description, copyleft, history and use added by Luca Montabone in June 2014.
;
;-----USE 
; Input values: u (3D zonal wind), v (3D meridional wind), ps (2D surface pressure), temp (3D temperature), thetapre (1D theta level array), lonar (1D longitude array in degrees, it could be a subset of total longitude array), latar (1D latitude array in degrees, it could be a subset of total latitude array), sigar (1D sigma level array), lonartot (1D all-longitude array in degrees), latartot (1D all-latitude array in degrees)
; Output values: PV (3D potential vorticity), utheta (3D zonal wind interpolated on theta levels), vtheta (3D meridional wind interpolated on theta levels), vort (relative vorticity on theta levels), dvdx, dudy (components of relative vorticity on theta levels).
; Note that it might be tempting to give the total / relative vorticity as an input to this script. However, this would be incorrect, as vorticity calculated on a particular level and subsequently interpolated onto theta surfaces is NOT necessarily the same as vorticity calculated by first interpolating the winds onto theta levels, and then calculating the gradients. 
;
; ------------------------------------- These variables are inputs---------------------------->|These variables are outputs-->|
PRO thetaint_version4_2015_earth, u,v,ps,temp,thetapre,lonar,latar,sigar,lonartot,latartot,PV,utheta,vtheta,vort,dvdx,dudy


; ------------------Set your planet's parameters----------------------------------


rm=6376.0e3 ;Planet radius in metres.
g=9.8 ;Graviational acceleration
Rct=287.04 ; Rgct is the R/g constant
Rgct=Rct/g
pref=1*1.0e5 ; Reference pressure for potential temperature (in Pascal).
r_cp=2./7.     ; R/cp value for power in potential temperature.
omega = 7.292e-5 ; Planet rotation rate omega, to go in the Coriolis parameter f, in units of radians per second. omega = 2*pi/rotation period.


; --------------------------------------------------------------------------------


; *****************************************************************************************************************************************************************************************************
; Description of criteria for longitude arrays and latitude arrays IF they are a SUBSET of the model's global coverage. (i.e. if the model goes from -87.5, 87.5 latitude, but you only pass data between -30,30, etc.)
; If subroutine gets passed a longitude array that is a subset of the true global dataset in longitude, the PV will not be calculated at the longitudes at the edge of the array, as it cannot calculate the relative vorticity using the central differencing method used in what follows. I.e. if you pass the routine an array of data in longitude corresponding to:
;
; lonar=[-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.]
;
; then PV will only be caluclated at longitudes:
; [-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.]
;
; So, if you wish to have PV calculated at +/-25 in this example, you'll need to provide data in the longitude range:
; lonar=[-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.]
;
;If, however, you pass data that is periodic in longitude, the PV will be calculated at all longitudes.
; 
; As the latitudinal direction is not periodic in the same way the longitudinal direction is, the PV will ALWAYS be outputted on an array that is smaller by 2 points in the latitudinal direction, even if the global dataset is passed. I.e. if
; latar=[-87.5,-85.,-82.5,......82.5,85.,87.5,]
; then the PV will only be outputted at:
; [-85.,-82.5,......82.5,85.]
;
; If a subset is passed in latitude, the same procedure is applied as is described for longitude above
; *****************************************************************************************************************************************************************************************************

;Need to see if prescribed theta levels are equally spaced.
dthetatest=(thetapre(n_elements(thetapre)-1)-thetapre(0))/(n_elements(thetapre)-1)
thetacomp=findgen(n_elements(thetapre))
thetacomp=thetapre(0)+(dthetatest*thetacomp)

if thetacomp(1) eq thetapre(1) then eqspace=1 else eqspace=0

if eqspace eq 0 then begin
print, 'theta levels not equally spaced'
stop
endif ; Nothing currently implemented to allow the use of non-equally spaced theta data.


;--------------------------------This is test to see if I am dealing with the entire global dataset, or a subset of it--------------------------------
 if latar(0) eq max(latartot) then correclatstart=0 else correclatstart=1 
 if latar(n_elements(latar)-1) eq min(latartot) then correclatend=0 else correclatend=1
 if n_elements(latar) eq n_elements(latartot) then wholelat=1 else wholelat=0

 if lonar(0) eq min(lonartot) then correclonstart=0 else correclonstart=1
 if lonar(n_elements(lonar)-1) eq max(lonartot) then correclonend=0 else correclonend=1
 if n_elements(lonar) eq n_elements(lonartot) then wholelon=1 else wholelon=0
;--------------------------------------------------------------------------------------------------------------------------------------------------------

;First thing is to interpolate the winds onto theta levels. To do this we have to calculate theta on the current model sigma levels...

;-------------------Initial theta calculation---------------------------------------------


prar=fltarr(n_elements(lonar),n_elements(latar),n_elements(sigar),n_elements(u(0,0,0,*)))

for k=0,n_elements(sigar)-1 do begin
  prar(*,*,k,*)=sigar(k)*ps(*,*,*) ;Simple calculation of 3D pressure. 
endfor

thetaar=temp*((pref)/prar)^(r_cp) ; Calculate theta on sigma levels, to enable interpolation onto theta levels.


print, 'interpolating winds'
utheta=fltarr(n_elements(lonar),n_elements(latar),n_elements(thetapre),n_elements(u(0,0,0,*)))
vtheta=fltarr(n_elements(lonar),n_elements(latar),n_elements(thetapre),n_elements(u(0,0,0,*)))
tempreal=fltarr(n_elements(lonar),n_elements(latar),n_elements(thetapre),n_elements(u(0,0,0,*)))

in_u=fltarr(n_elements(thetapre))
in_v=fltarr(n_elements(thetapre))
in_t=fltarr(n_elements(thetapre))
thetalev=fltarr(n_elements(thetapre))


  for nt=0,(n_elements(u(0,0,0,*))-1) do begin
    for nlat=0,(n_elements(u(0,*,0,0))-1) do begin
      for nlon=0,(n_elements(u(*,0,0,0))-1) do begin

	psurf=ps(nlon,nlat,nt)
	tprof=temp(nlon,nlat,*,nt)
	uprof=u(nlon,nlat,*,nt)
	vprof=v(nlon,nlat,*,nt)
	thetalev=thetaar(nlon,nlat,*,nt)



	  for nlev=0L,(n_elements(thetapre)-1) do begin
	  WEI_VERT,0,thetalev,thetapre(nlev),lolev,vw ; need to do vertical weight calculation once for interpolating the data on thetalev levels onto thetapre levels

	    if(lolev ge 0) then begin

	    in_t(nlev)=tprof(lolev)+(tprof(lolev+1)-tprof(lolev))*vw
	    in_u(nlev)=uprof(lolev)+(uprof(lolev+1)-uprof(lolev))*vw
	    in_v(nlev)=vprof(lolev)+(vprof(lolev+1)-vprof(lolev))*vw

	    endif else begin
	    in_t(nlev)=0
	    in_u(nlev)=0
	    in_v(nlev)=0

	  endelse

  utheta(nlon,nlat,*,nt)=in_u(*)
  tempreal(nlon,nlat,*,nt)=in_t(*)
  vtheta(nlon,nlat,*,nt)=in_v(*)


	  endfor

	endfor
      endfor
;     print, nt
    endfor

;--------------------------------------Now have winds on Theta Levels-----------------------------------------------------
;---------------------------- Now want to calculate Vorticity on Theta Levels---------------------------------------------

;This bit does the absolute vorticity calculations.
print, 'Calculating Absolute Vorticity'
dvdx=fltarr(n_elements(u(*,0,0,0)),n_elements(u(0,*,0,0)),n_elements(thetapre),n_elements(u(0,0,0,*)))
dudy=fltarr(n_elements(u(*,0,0,0)),n_elements(u(0,*,0,0)),n_elements(thetapre),n_elements(u(0,0,0,*)))
vort=fltarr(n_elements(u(*,0,0,0)),n_elements(u(0,*,0,0)),n_elements(thetapre),n_elements(u(0,0,0,*)))

pi= 3.141592653589793

latarpi=((2*pi)/(360.0))*latar
lonarpi=((2*pi)/(360.0))*lonar

; X-Gradient **********************************************************************************
print,'running x grad'
for lat=correclatstart,(n_elements(latar)-1-correclatend) do begin

 for i=1,(n_elements(lonar)-1-1) do begin ; This central-difference method works out gradient at point a using point a+1 and point a-1, not point a itself! :)
   dvdx(i,lat,*,*)=(vtheta(i+1,lat,*,*)-vtheta(i-1,lat,*,*))/(((lonarpi(i+1)*rm*cos(latarpi(lat)))-(lonarpi(i-1)*rm*cos(latarpi(lat))))) 
 endfor

 if wholelon then begin ; Can only calculate the dvdx at the longitudes at the edge of the domain if it is periodic. wholelon=1 if longitudes are periodic.
   ; x-gradient for first longitude (0)
   dvdx(0,lat,*,*)=(vtheta(1,lat,*,*)-vtheta((n_elements(lonar)-1),lat,*,*))/(((lonarpi(2)-lonarpi(0))*rm*cos(latarpi(lat)))) ; uses lonarpi(2)-lonarpi(0) as dlambda is constant anyway, but also that lonarpi(1)-lonarpi(n_elements(lonar)-1) would be ~360. bc of minus signs. 

   ; x-gradient for last longitude (n_elements(lonar)-1)
   dvdx((n_elements(lonar)-1),lat,*,*)=(vtheta(0,lat,*,*)-vtheta((n_elements(lonar)-2),lat,*,*))/(((lonarpi(2)-lonarpi(0))*rm*cos(latarpi(lat))))  
 endif
;print, lat
endfor

; Y-Gradient **************************************************************************************************
print,'running y grad'

for lon=correclonstart,(n_elements(lonar)-1-correclonend) do begin 

 for i=1,(n_elements(latar)-1-1) do begin
   dudy(lon,i,*,*)=(utheta(lon,i+1,*,*)-utheta(lon,i-1,*,*))/(((latarpi(i+1)*rm)-(latarpi(i-1)*rm)))-(utheta(lon,i,*,*)*tan(latarpi(i))/rm) ; For origin of second term see Middle Atmosphere Dynamics (Andrews et al) 1987 equation 3.8.4b. Minus sign to make sure that when we do dv/dx-du/dy we have the correct sign for the + u*tan(lat)/rm 
 endfor
endfor

vort=dvdx-dudy ; At this point, the relative vorticity is defined everywhere on potential temperature levels.


; ******************************************************Begin enforced sub-setting of y data**************************************************************
; At this point, because we have not calculated the y gradient of u at the two edge latitudinal gridpoints, and therefore we cannot claim to have a vorticity, or a PV defined there. So, we throw all of the data at the two final gridpoints in latitude away. Do that with the whole dataset, so as to make the rest of the routine think it is simply dealing with a dataset of nlat-2 gridpoints in y, as opposed to one with nlat.
; stop


vort=vort(*,1:(n_elements(vort(0,*,0,0))-2),*,*)
dvdx=dvdx(*,1:(n_elements(dvdx(0,*,0,0))-2),*,*)
dudy=dudy(*,1:(n_elements(dudy(0,*,0,0))-2),*,*)
utheta=utheta(*,1:(n_elements(utheta(0,*,0,0))-2),*,*)
vtheta=vtheta(*,1:(n_elements(vtheta(0,*,0,0))-2),*,*)
v=v(*,1:(n_elements(v(0,*,0,0))-2),*,*)
u=u(*,1:(n_elements(u(0,*,0,0))-2),*,*)
prar=prar(*,1:(n_elements(prar(0,*,0,0))-2),*,*)
ps=ps(*,1:(n_elements(ps(0,*,0,0))-2),*)
temp=temp(*,1:(n_elements(temp(0,*,0,0))-2),*,*)
thetaar=thetaar(*,1:(n_elements(thetaar(0,*,0,0))-2),*,*)
tempreal=tempreal(*,1:(n_elements(tempreal(0,*,0,0))-2),*,*)
latar=latar(1:(n_elements(latar)-2))
latarpi=latarpi(1:(n_elements(latarpi)-2))

; stop

; ******************************************************End enforced sub-setting of y data**************************************************************


; ******************************************************Begin optional sub-setting of x data**************************************************************
; If we are given the whole longitudinal extent of the data (i.e. data that covers the whole globe in the longitudinal direction), then x gradients at the two 'edges' of the domain (i.e. lonar(0) and lonar(-1)), can be calculated using this periodicity in x. So, there is no enforced sub-setting of the data in the x direction, as there is in y. However, if the whole longitudinal dataset is not passed to the routine (i.e. wholelon=0) then the x gradient of the data cannot be calculated at the extrema in x, and so the data must be sub-sampled in x.

if ~wholelon then begin
vort=vort(1:(n_elements(vort(*,0,0,0))-2),*,*,*)
dvdx=dvdx(1:(n_elements(dvdx(*,0,0,0))-2),*,*,*)
dudy=dudy(1:(n_elements(dudy(*,0,0,0))-2),*,*,*)
utheta=utheta(1:(n_elements(utheta(*,0,0,0))-2),*,*,*)
vtheta=vtheta(1:(n_elements(utheta(*,0,0,0))-2),*,*,*)
u=u(1:(n_elements(u(*,0,0,0))-2),*,*,*)
v=v(1:(n_elements(v(*,0,0,0))-2),*,*,*)
prar=prar(1:(n_elements(prar(*,0,0,0))-2),*,*,*)
ps=ps(1:(n_elements(vort(*,0,0,0))-2),*,*)
temp=temp(1:(n_elements(vort(*,0,0,0))-2),*,*,*)
thetaar=thetaar(1:(n_elements(thetaar(*,0,0,0))-2),*,*,*)
tempreal=tempreal(1:(n_elements(tempreal(*,0,0,0))-2),*,*,*)
lonar=lonar(1:(n_elements(lonar)-2))
lonarpi=lonarpi(1:(n_elements(lonarpi)-2))
endif
; ******************************************************End optional sub-setting of x data**************************************************************



f=2*(omega)*sin(latarpi)
vorttot=fltarr(n_elements(vort(*,0,0,0)),n_elements(vort(0,*,0,0)),n_elements(vort(0,0,*,0)),n_elements(vort(0,0,0,*)))

for ii=0,(n_elements(latar)-1) do begin
 vorttot(*,ii,*,*)=vort(*,ii,*,*)+f(ii)
endfor

;Up to now I have calculated total vorticity on theta levels and now we need dtheta/dpressure on those same theta levels. For us to do this simply would require that the theta levels are equally spaced.



print, 'Now for the potential temperature bit.......'

;What happens here is that the dtheta/dpressure part of the PV equation is worked out. It is done by calculating, point by point, through the entire dataset, dtheta/dp based on the data grided on sigma levels. But we know that because the sigma levels are not regularly gridded in the vertical, the derrivatives that we are forced to calculate are central difference derrivatives, but these derrivatives are calculated on new sigma levels in-between the current sigma levels.

; So say we have 2 sigma levels, sig0,sig1 :



;_______________________________________________________ sig1


;-------------------------------------------------------- sig0.5

;________________________________________________________ sig0


;To work out dtheta/dp imagine we make a new sigma level where sig0.5=(sig1+sig0)*0.5 Then dtheta/dp(sigma=sig0.5)=(theta(sig1)-theta(sig0))/(2.0*(pr(sig1)-pr(sig0))) 
; This is the only way we can get sigar to allow us to calculate dtheta/dp using centrally differenced derrivatives. Once this is done, we have a dtheta/dp array on these new sigma levels. The new sigma levels have a thetaarr associated with them, which is then used to interpolate the results onto thetapre.


dthetaarr=fltarr(n_elements(vorttot(*,0,0,0)),n_elements(vorttot(0,*,0,0)),(n_elements(sigar)-1),n_elements(vorttot(0,0,0,*)))
sigarnew=fltarr(n_elements(sigar)-1)
dpresarr=fltarr(n_elements(vorttot(*,0,0,0)),n_elements(vorttot(0,*,0,0)),(n_elements(sigar)-1),n_elements(vorttot(0,0,0,*)))
temptwo=fltarr(n_elements(vorttot(*,0,0,0)),n_elements(vorttot(0,*,0,0)),(n_elements(sigar)-1),n_elements(vorttot(0,0,0,*)))
prartwo=fltarr(n_elements(vorttot(*,0,0,0)),n_elements(vorttot(0,*,0,0)),(n_elements(sigar)-1),n_elements(vorttot(0,0,0,*)))


for jj=1,(n_elements(sigar)-1) do begin

; sigarnew(jj-1)=exp((alog(sigar(jj))+alog(sigar(jj-1))/2.0)) ; This would be one way of calculating the new sigma level, but it's not what's required, as we are not looking for the level that is midway between in terms of height, but in terms of sigma. Hence we use the below instead:

sigarnew(jj-1)=((sigar(jj))+(sigar(jj-1)))/2.0 ; Works out sigma of new levels!

dpresarr(*,*,jj-1,*)=(prar(*,*,jj,*)-prar(*,*,jj-1,*))
dthetaarr(*,*,jj-1,*)=(thetaar(*,*,jj,*)-thetaar(*,*,jj-1,*))/(dpresarr(*,*,jj-1,*))
temptwo(*,*,jj-1,*)=(temp(*,*,jj,*)+temp(*,*,jj-1,*))/2.0 ; Works out the temperature on the new levels using linear interpolation. Is simple because we are just looking for the level halfway between the two existing levels.

endfor

for pr=0, (n_elements(sigarnew)-1) do begin
prartwo(*,*,pr,*)=sigarnew(pr)*ps ; Calculate new pressure levels based on new sigma levels.
endfor


thetaartwo=temptwo*((pref)/prartwo)^(r_cp)
thetalevtwo=fltarr(n_elements(thetaartwo(0,0,*,0)))
in_var=fltarr(n_elements(thetapre))


dthetaarrth=fltarr(n_elements(vorttot(*,0,0,0)),n_elements(vorttot(0,*,0,0)),(n_elements(thetapre)),n_elements(vorttot(0,0,0,*)))

  for nt=0,(n_elements(dthetaarr(0,0,0,*))-1) do begin
    for nlat=0,(n_elements(dthetaarr(0,*,0,0))-1) do begin
      for nlon=0,(n_elements(dthetaarr(*,0,0,0))-1) do begin

	varprof=dthetaarr(nlon,nlat,*,nt)
 	thetalevtwo=thetaartwo(nlon,nlat,*,nt)
	  for nlev=0L,(n_elements(thetapre)-1) do begin
	  WEI_VERT,0,thetalevtwo,thetapre(nlev),lolevtwo,vwtwo ; need to do vertical weight calculation once for interpolating the old thetaarr onto new theta levels
	    if(lolevtwo ge 0) then begin

	    in_var(nlev)=varprof(lolevtwo)+(varprof(lolevtwo+1)-varprof(lolevtwo))*vwtwo
;

	    endif else begin
	    in_var(nlev)=0
	    endelse


  dthetaarrth(nlon,nlat,*,nt)=in_var(*)


 	  endfor

	endfor
      endfor
     print, nt
    endfor

;Final PV values
PV=-g*vorttot*dthetaarrth


print, 'saving'
save,filename=filenamestr
end
