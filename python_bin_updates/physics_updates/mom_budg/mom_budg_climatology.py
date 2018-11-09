"""
Evaluate and plot momentum budget at 150 hPa - redo in line with reviewer 1s suggestions, see how it looks

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling_updates import month_dic, gradients as gr
import sh
from pylab import rcParams
from mom_budg import mom_budg_fun
    
    
    
def mom_budg_clim(run, lev=150, plot_list=['fv_ageo'], land_mask=None):
    
    data = mom_budg_fun(run, lev=lev)
    # Take seasonal averages
    data.coords['season'] = np.mod(data.xofyear + 5., 72.) // 18. 
    #print data.season
    data = data.groupby('season').mean(('xofyear')).sel(pfull=lev)
    
    rcParams['figure.figsize'] = 10, 8
    rcParams['font.size'] = 16
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/mom_budg/climatology/' + run + '/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    levels = np.arange(-20,21.1,2.)
    
    for var in plot_list:
        
        # Four subplots
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = data[var][i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = levels, add_labels=False, add_colorbar=False, extend='both')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask==None:
                land = xr.open_dataset(land_mask)
                land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        plt.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.01, hspace=0.3, wspace=0.3)
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.1, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + var + '_' + run + '.pdf', format='pdf')
        plt.close()



def mom_budg_monthly(run, lev=150, plot_list=['fv_ageo'], land_mask=None):
    
    data = mom_budg_fun(run, lev=lev)
    # Take month averages
    data.coords['month'] = (data.xofyear -1) // 6 
    data = data.groupby('month').mean(('xofyear')).sel(pfull=lev)
    
    rcParams['figure.figsize'] = 10, 8
    rcParams['font.size'] = 16
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/mom_budg/climatology/' + run + '/monthly/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
        
    levels = np.arange(-20,21.1,2.)
    
    for var in plot_list:
        
        # Four subplots
        f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
        title = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug','Sep', 'Oct', 'Nov', 'Dec']
        
        
        for i in range(12):            
            f1 = data[var][i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = levels, add_labels=False, add_colorbar=False, extend='both')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask==None:
                land = xr.open_dataset(land_mask)
                land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
            axes[i].set_xticks(range(0,361,120))
            axes[i].set_yticks(range(-60,61,30))
            axes[i].set_ylim([-60,60])
        
        for ax in [ax1,ax5,ax9]:
            ax.set_ylabel('Latitude')
        for ax in [ax9,ax10,ax11,ax12]:
            ax.set_xlabel('Longitude')
        
        plt.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.01, hspace=0.3, wspace=0.3)
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.1, pad=0.15, aspect=30, shrink=0.5)
                
        plt.savefig(plot_dir + var + '_' + run + '.pdf', format='pdf')
        plt.close()
        
        

if __name__ == "__main__":
    
    plot_list = ['mom_cross', 'mom_trans', 'mom_stat', 'fv_local', 'fv_ageo', 'dphidx', 'mom_sum']
    mom_budg_clim('half_nh_shallow', plot_list=plot_list, land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_shallow.nc')
    mom_budg_monthly('half_nh_shallow', plot_list=plot_list, land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_shallow.nc')
    
    mom_budg_clim('control_qflux', plot_list=plot_list, land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    mom_budg_monthly('control_qflux', plot_list=plot_list, land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    
    