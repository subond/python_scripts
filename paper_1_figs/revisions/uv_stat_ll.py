"""
Evaluate and plot u*v* for full expt for JJA at 150 hPa cf fig 5 Shaw 2014 

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from pylab import rcParams

    
    
def uv_stat_ll(run):
    
    rcParams['figure.figsize'] = 6, 8
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/paper_1_figs/revisions/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/'+run+'.nc')
    
    u_jja = data.ucomp.sel(xofyear=range(31,49)).mean('xofyear')
    v_jja = data.vcomp.sel(xofyear=range(31,49)).mean('xofyear')
    
    #uv_stat_jja = u_jja * v_jja - u_jja.mean('lon') * v_jja.mean('lon')
    
    uv_stat_jja_mima = (u_jja - u_jja.mean('lon')) * (v_jja - v_jja.mean('lon'))
    
    #uv_stat_jja.sel(pfull=150.).plot.contourf(levels=np.arange(-200.,201.,20.))
    

    #plt.show()
    
    data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/era_mom_vars.nc')
    
    u_jja = data.ucomp.sel(day_of_yr=range(152,244)).mean('day_of_yr')
    v_jja = data.vcomp.sel(day_of_yr=range(152,244)).mean('day_of_yr')
    
    uv_stat_jja_era = (u_jja - u_jja.mean('lon')) * (v_jja - v_jja.mean('lon'))
    
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    
    uv_stat_jja_era.sel(pfull=150.).plot.contourf(x='lon', y='lat', ax=ax1, levels=np.arange(-150.,151.,20.), extend='both', add_labels=False, add_colorbar=False)
    land.land_mask.plot.contour(x='lon', y='lat', ax=ax1, levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k')
    ax1.set_ylabel('Latitude')
    ax1.set_ylim(-60,60)
    ax1.set_yticks(np.arange(-60.,61.,30.))
    ax1.set_xticks(np.arange(0.,361.,60.))
    ax1.grid(True,linestyle=':')
    ax1.text(-60, 60, 'a)')
    
    f1 = uv_stat_jja_mima.sel(pfull=150.).plot.contourf(x='lon', y='lat', ax=ax2, levels=np.arange(-150.,151.,20.), extend='both', add_labels=False, add_colorbar=False)
    land.land_mask.plot.contour(x='lon', y='lat', ax=ax2, levels = np.arange(-5.,7.,5.), add_colorbar=False, add_labels=False, colors='k')
    ax2.set_ylabel('Latitude')
    ax2.set_xlabel('Longitude')
    ax2.set_ylim(-60,60)
    ax2.set_yticks(np.arange(-60.,61.,30.))
    ax2.set_xticks(np.arange(0.,361.,60.))
    ax2.grid(True,linestyle=':')
    ax2.text(-60, 60, 'b)')
    
    
    plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
    #Colorbar
    cb1=fig.colorbar(f1, ax=(ax1, ax2), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.12, aspect=30, ticks=np.arange(-150,151,60))
    #cb1.set_label('Stationary eddy momentum transport, m$^2$s$^{-2}$')
 
    plt.savefig(plot_dir+'uv_stat_ll.pdf', format='pdf')
    plt.close()
    


uv_stat_ll('full_qflux')



