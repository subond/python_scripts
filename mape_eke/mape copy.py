# Load up climatology data and calculate MAPE

import xarray as xr
from physics import model_constants as mc, gradients as gr
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as spint

def trop_height(temp):
    # Use temperature to calculate tropopause level
    # Take vertical gradient wrt height and find pressure level at which dT/dz goes and stays above -2 K/km
    
    levs = temp.pfull[temp.pfull <=500.]
    
    pfull_big = np.arange(1000.,0.,-10.)
    
    f = spint.interp1d(temp.pfull, temp.mean('lon'), axis=1, fill_value='extrapolate')
    temp_big = f(pfull_big)
    temp_big = xr.DataArray(temp_big, coords=[temp.xofyear, pfull_big, temp.lat], dims=['xofyear', 'pfull', 'lat'])
    
    rho_big = temp_big.pfull * 100. / mc.rdgas / temp_big
    dtdp_big = gr.ddp(temp_big)
    levs_big = dtdp_big.pfull[dtdp_big.pfull <=500.]
    dtdz_big = (-1. * dtdp_big * rho_big * mc.grav * 1000.).sel(pfull=levs_big)
        
    
    rho = temp.pfull * 100. / mc.rdgas / temp.mean('lon')
    dtdp = gr.ddp(temp.mean('lon'))
    dtdz = (-1. * dtdp * rho * mc.grav * 1000.)#.sel(pfull=levs)
    
    f = spint.interp1d(dtdz.pfull, dtdz,  axis=1, fill_value='extrapolate')
    dtdz_big = f(levs_big)
    
        
   # dtdz_mask = xr.DataArray( np.where(dtdz >= -2., 1., 0.), coords=[dtdz.xofyear, dtdz.pfull, dtdz.lat], dims=['xofyear', 'pfull', 'lat'])
    dtdz_mask = np.where(dtdz >= -2., 1., 0.)
    
    #print dtdz_mask[0,:-1,0]
    #print dtdz_mask[0,1:,0]
    
    #print (dtdz_mask[0,:-1,0] != dtdz_mask[0,1:,0])
    
    #print np.where(dtdz_mask[0,:-1,0] != dtdz_mask[0,1:,0])
    
    trop = np.zeros([len(dtdz.xofyear), len(dtdz.lat)])
    for i in range(len(dtdz.xofyear)):
        for j in range(len(dtdz.lat)):
            #trop[i,j] = min(dtdz.pfull[np.where(dtdz[i,:,j] <= -2.)])
            trop[i,j] = min(levs_big[np.where(dtdz_big[i,:,j] <= -2.)])
#            trop[i,j] = dtdz.pfull[ np.where(dtdz_mask[i,:-1,j] != dtdz_mask[i,1:,j])[0] ]
    
    plt.contourf(trop, levels=np.arange(90.,260.,10.))
    plt.colorbar()
    
    print trop[20,:]
    print dtdz[0,:,0]
    print np.where(dtdz[0,:,0] >= -2.)
    
    plt.figure(2)
    plt.contourf(dtdz_big[18,:,:], levels=np.arange(-12.,13.,1.))
    plt.show()
    
    
    
  #  for i=2:89
  #  trop_p_000(i) = rC(min(find(dtdz_000_zav(i,:) >= -2) )) ;
  #  end
    
    
    

def mape(run):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    theta = data.temp * (1000./data.pfull)**mc.kappa
    
    trop_height(data.temp)
    
    
    # Evaluate tropopause height - create new function for this



#thetall=cube2latlon(xc,yc,theta_tav,xi,yi);
#theta_zav(:,:) = mean(thetall,1);
#thetasq_zav = theta_zav.^2;
#dthetadp = gradient(theta_zav,-4000);
#for i=1:90
#dthetadp_weight(i,:) = cos(yi(i).*pi./180) .* dthetadp(i,:) ;
#theta_weight(i,:) = cos(yi(i).*pi./180) .* theta_zav(i,:) ;
#thetasq_weight(i,:) = cos(yi(i).*pi./180) .* thetasq_zav(i,:) ;
#end

#dthetadp_mean_n = sum(dthetadp_weight(x000n-7:x000n+7,:),1)./sum(cos(yi(x000n-7:x000n+7).*pi./180));
#dthetadp_mean_s = sum(dthetadp_weight(x000s-7:x000s+7,:),1)./sum(cos(yi(x000s-7:x000s+7).*pi./180));
#dthetadp_mean = (dthetadp_mean_n+dthetadp_mean_s)./2;

#theta_mean_n = sum(theta_weight(x000n-7:x000n+7,:),1)./sum(cos(yi(x000n-7:x000n+7).*pi./180));
#theta_mean_s = sum(theta_weight(x000s-7:x000s+7,:),1)./sum(cos(yi(x000s-7:x000s+7).*pi./180));
#theta_mean = (theta_mean_n + theta_mean_s)./2;

#thetasq_mean_n = sum(thetasq_weight(x000n-7:x000n+7,:),1)./sum(cos(yi(x000n-7:x000n+7).*pi./180));
#thetasq_mean_s = sum(thetasq_weight(x000s-7:x000s+7,:),1)./sum(cos(yi(x000s-7:x000s+7).*pi./180));
#thetasq_mean = (thetasq_mean_n + thetasq_mean_s)./2;

#for k=1:25
#gamma_000(k) = -kappa./rC(k) .* 1./dthetadp_mean(k);
#end
#theta_var_000 = thetasq_mean - theta_mean.^2;
#mape_integrand_000 = (rC'./p0).^kappa .*gamma_000 .* theta_var_000 ;
#mape_000 = prefac .* sum(mape_integrand_000(3:trop_lev_000(x000n))) .* 4000./p0;
#%mape_000 = prefac .* sum(mape_integrand_000(b_lev_000(x000n) :trop_lev_000(x000n))) .* 4000./p0./(trop_lev_000(x000n)- b_lev_000(x000n) +1);

mape('full_qflux')
