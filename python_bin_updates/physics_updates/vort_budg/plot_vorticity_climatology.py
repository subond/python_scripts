"""
14/02/2018
Plot vortex stretching, horizontal advection, and absolute vorticity for a given run, or their differences for 2 runs

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from pylab import rcParams
from plot_vorticity_breakdown import vort_budg_terms


plot_dir = '/scratch/rg419/plots/vorticity_eq_clean/climatology/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

def vort_climatology(run, lev=150., stretching=True, horiz_adv=True, eddies=True, vor=True, land_mask=None):
    '''Plot vorticity budget ll. RG 
       Imputs: run = run_name
               lev = level to plot
               lonin = longitude range to average over
               planetary_only = only plot planetary vorticity terms
               month_labels = label x axis with month of year (no good for non earth year length)
               rot_fac = scale factor for Earth's rotation rate
               no_eddies = don't plot transient eddies too'''
    
    rcParams['figure.figsize'] = 10, 8
    rcParams['font.size'] = 16
    
    
    # Call the above function to get fields to plot. For a Hovmoller, run is by definition not steady state. Pass other inputs on.
    ds = vort_budg_terms(run, ll=True)
    
    # Take seasonal averages
    ds.coords['season'] = np.mod(ds.xofyear + 5., 72.) // 18. 
    #print data.season
    ds = ds.groupby('season').mean(('xofyear')).sel(pfull=lev)
    
    plt.set_cmap('RdBu_r')
    
    levels = np.arange(-1.5,1.6,0.25)
    
    
    if stretching:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = ds.stretching[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = levels, add_labels=False, add_colorbar=False, extend='both')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask==None:
                land = xr.open_dataset(land_mask)
                land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'stretching_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    if horiz_adv:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = ds.horiz_md[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = levels, add_labels=False, add_colorbar=False, extend='both')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask==None:
                land = xr.open_dataset(land_mask)
                land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'horiz_adv_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    
    if eddies:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = ds.transient_hm[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = levels, add_labels=False, add_colorbar=False, extend='both')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask==None:
                land = xr.open_dataset(land_mask)
                land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'eddies_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    
    if vor:
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        axes = [ax1,ax2,ax3,ax4]
        title = ['DJF', 'MAM', 'JJA', 'SON']
        
        for i in range(4):            
            f1 = ds.vor[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels=np.arange(-12.,13.,2.), add_labels=False, add_colorbar=False, extend='both')
            axes[i].grid(True,linestyle=':')
            axes[i].set_title(title[i])
            if not land_mask==None:
                land = xr.open_dataset(land_mask)
                land.land_mask.plot.contour(x='lon', y='lat', ax=axes[i], levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        
        for ax in [ax1,ax3]:
            ax.set_ylabel('Latitude')
        for ax in [ax3,ax4]:
            ax.set_xlabel('Longitude')
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'abs_vort_' + run + '.pdf', format='pdf')
        plt.close()






def vort_climatology_monthly(run, lev=150., stretching=True, horiz_adv=True, eddies=True, vor=True, land_mask=None):
    '''Plot vorticity budget Hovmoller. RG 3/11/2017
       Imputs: run = run_name
               lev = level to plot
               lonin = longitude range to average over
               planetary_only = only plot planetary vorticity terms
               month_labels = label x axis with month of year (no good for non earth year length)
               rot_fac = scale factor for Earth's rotation rate
               no_eddies = don't plot transient eddies too'''
    
    rcParams['figure.figsize'] = 10, 8
    rcParams['font.size'] = 16
    
    
    # Call the above function to get fields to plot. For a Hovmoller, run is by definition not steady state. Pass other inputs on.
    ds = vort_budg_terms(run, ll=True)
    
    # Take month averages
    ds.coords['month'] = (ds.xofyear -1) // 6 
    ds = ds.groupby('month').mean(('xofyear')).sel(pfull=lev)
    
    plt.set_cmap('RdBu_r')
    
    levels = np.arange(-1.5,1.6,0.25)
    
    title = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug','Sep', 'Oct', 'Nov', 'Dec']
    
    if stretching:
        f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
        
        for i in range(12):            
            f1 = ds.stretching[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = levels, add_labels=False, add_colorbar=False, extend='both')
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
        
        plt.savefig(plot_dir + 'stretching_monthly_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    if horiz_adv:
        f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
        
        for i in range(12):            
            f1 = ds.horiz_md[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = levels, add_labels=False, add_colorbar=False, extend='both')
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
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'horiz_adv_monthly_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    
    if eddies:
        f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
        
        for i in range(12):            
            f1 = ds.transient_hm[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels = levels, add_labels=False, add_colorbar=False, extend='both')
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
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'eddies_monthly_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    
    if vor:
        f, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12)) = plt.subplots(3, 4, sharex='col', sharey='row')
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12]
        
        for i in range(12):            
            f1 = ds.vor[i,:,:].plot.contourf(x='lon', y='lat', ax=axes[i], levels=np.arange(-12.,13.,2.), add_labels=False, add_colorbar=False, extend='both')
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
        
        cb1=f.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.15, aspect=30, shrink=0.5)
        
        plt.savefig(plot_dir + 'abs_vort_monthly_' + run + '.pdf', format='pdf')
        plt.close()
    
    
    

    
if __name__ == "__main__":
        
    #vort_climatology('idealised_2cont', land_mask = '/scratch/rg419/Experiments/wavenumber_2/input/1hill.nc')
    #vort_climatology('idealised_2am')
    #vort_climatology('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')    
    vort_climatology_monthly('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')    
    vort_climatology_monthly('half_nh_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_shallow.nc')    
    #vort_climatology('control_qflux', land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')    
    #vort_climatology('no_americas', land_mask = '/scratch/rg419/python_scripts/land_era/land_era_no_america.nc')
    #vort_climatology('frozen_am', land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    #vort_climatology('no_TIP', land_mask = '/scratch/rg419/Isca/input/land_masks/era_land_t42.nc')
    
    
    
    
    
    