import numpy as np

def wei_vert(linlog,vertlev,intvert):

	lev=len(vertlev)

	lolev=-1
	vw=0

	k_tick=range(lev-1)

	for k in k_tick:
		if (vertlev[k+1]-vertlev[k]) > 0.0:
			a=vertlev[k]
			b=vertlev[k+1]
		else:
			a=vertlev[k+1]
			b=vertlev[k]
		if intvert >= a and intvert <= b:
			lolev=k
			if linlog ==0:
				vw=np.log(intvert/vertlev[k])/np.log(vertlev[k+1]/vertlev[k])			
			elif linlog==1:
				vw=(intvert - vertlev[k])/(vertlev[k+1] - vertlev[k])			
	
	
	return lolev,vw

#PRO wei_vert,linlog,vertlev,intvert,lolev,vw
#
#; WEI_VERT calculates the log weight for pressure coordinates (0) or the
#; linear weight for altitude coordinates (1)
#
#lev=n_elements(vertlev)
#
#lolev=-1
#for k=0,lev-2 do begin 
#  if((vertlev(k+1)-vertlev(k)) gt 0.0) then begin
#    a=vertlev(k)
#    b=vertlev(k+1)
#  endif else begin 
#    a=vertlev(k+1)
#    b=vertlev(k)
#  endelse
#  if(intvert ge a and intvert le b) then begin
#    lolev=k
#    if(linlog eq 0) then vw=alog(intvert/vertlev(k))/alog(vertlev(k+1)/vertlev(k))
#    if(linlog eq 1) then vw=(intvert-vertlev(k))/(vertlev(k+1)-vertlev(k))    
#  endif
#endfor

#end
