# Compare dv/dy with vmax/half cell width

"""
Load in the vorticity budget terms and produce a lat-time plot at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import month_dic, cell_area
import sh
from physics import gradients as gr
from pylab import rcParams

def rad_eq_t(run, pentad, lev=150, period_fac=1.):
    
    rcParams['figure.figsize'] = 15, 6.25
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/crit_lat_test/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    area = cell_area(42, '/scratch/rg419/GFDL_model/GFDLmoistModel/')
    
    #Load in data
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    data['area'] = (('lat','lon'), area)
    
    # Radiative equilibrium temperature
    stefan = 5.6734e-8
    t_rad_eq = (data.toa_sw/stefan) ** (1./4.)
    
    mn_dic = month_dic(1)
    tickspace = np.arange(13,72,18) * period_fac
    labels = [mn_dic[(k+5)/6 ] for k in range(13, 72, 18)]
    levels = np.arange(-1.5,1.6,0.25)

    t_rad_eq.mean('lon').plot.contourf(x='xofyear', y='lat', extend = 'both', add_labels=False, levels=np.arange(0.,301.,10.))
    plt.ylabel('Latitude')
    plt.xlabel('')
    plt.yticks(np.arange(-60,61,30))
    plt.xticks(tickspace,labels,rotation=25)
    plt.title('T, K', fontsize=17)
    plt.grid(True,linestyle=':')
    plt.tight_layout()  
    
    figname = 't_rad_eq_' + run + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    t_rad_eq.mean(('lon', 'xofyear')).plot()
    plt.ylabel('Latitude')
    plt.xlabel('Radiative equilibrium temperature')
    figname = 't_rad_eq_mean_' + run + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    data.t_surf.mean('lon').plot.contourf(x='xofyear', y='lat', extend = 'both', add_labels=False, levels=np.arange(240.,311.,5.))
    plt.ylabel('Latitude')
    plt.xlabel('')
    plt.yticks(np.arange(-60,61,30))
    plt.xticks(tickspace,labels,rotation=25)
    plt.title('T, K', fontsize=17)
    plt.grid(True,linestyle=':')
    plt.tight_layout()  
    figname = 't_surf_' + run + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    delta_t_N = data.t_surf.mean('lon').max('lat') - data.t_surf.mean('lon').isel(lat=range(32,64)).min('lat')
    delta_t_S = data.t_surf.mean('lon').max('lat') - data.t_surf.mean('lon').isel(lat=range(0,32)).min('lat')
    t_mean = ((data.t_surf * data.area).sum(('lat','lon'))/data.area.sum(('lat','lon')) ).mean('xofyear')
    
    print run, t_mean.values, delta_t_S[pentad].values
    
    delta_t_N.plot(color='k')
    delta_t_S.plot(color='b')
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.ylabel('delta T, K')
    plt.ylim([5,65])
    plt.grid(True,linestyle=':')
    plt.tight_layout()  
    figname = 'deltaT_' + run + '.pdf'
    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
    

rad_eq_t('sn_1.000', pentad=46)
rad_eq_t('sn_2.000', pentad=85, period_fac=2.)
rad_eq_t('sn_0.500', pentad=26, period_fac=0.5)
rad_eq_t('rt_2.000', pentad=59)
rad_eq_t('rt_0.500', pentad=46)