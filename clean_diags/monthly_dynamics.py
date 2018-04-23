"""
Functions to plot
1) Streamfunction and velocity potential
2) Kinetic energy associated with irrotational and nondivergent motions
3) Absolute vorticity

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, cell_area, month_dic
from windspharm.xarray import VectorWind
import sh

def sf_vp(run, months, filename='plev_pentad', timeav='month', period_fac=1.,land_mask=False,level=9):
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/sf_vp/'+str(level)+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    mn_dic = month_dic(1)
    
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    uwnd = data.ucomp[:,level,:,:]
    vwnd = data.vcomp[:,level,:,:]

    w = VectorWind(uwnd, vwnd)
    sf, vp = w.sfvp()
    sf *= 1e-6
    vp *= 1e-6
    
    # Plot streamfunction.
    for i in range(0,12):
        ax=sf[i,:,:].plot.contourf(x='lon',y='lat',levels=np.arange(-180.,181.,20.), extend='both', add_colorbar=False, add_label = False)
        cb1=plt.colorbar(ax)
        cb1.set_label('Streamfunction')
        # Plot velocity potential.
        if land_mask:
            land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')
            land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
            land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='0.75',add_colorbar=False,add_labels=False)
        cs = vp[i,:,:].plot.contour(x='lon',y='lat', levels=np.arange(-16.,16.1,4.), extend='both', add_colorbar=False, colors='k', add_label = False)    
        plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.ylim(-10,45)
        plt.xlim(25,150)
        plt.tight_layout()  
        plt.savefig(plot_dir+str(i+1)+'_'+str(mn_dic[i+1])+'.png')
        plt.close()



def ke_components(run, months, filename='plev_pentad', timeav='month', period_fac=1.,land_mask=False,level=9):
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/ke/'+str(level)+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    mn_dic = month_dic(1)
    
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    uwnd = data.ucomp[:,level,:,:]
    vwnd = data.vcomp[:,level,:,:]

    w = VectorWind(uwnd, vwnd)
    uchi, vchi, upsi, vpsi = w.helmholtz()
    ke_chi = 0.5*(uchi*uchi + vchi*vchi)
    ke_psi = 0.5*(upsi*upsi + vpsi*vpsi)

    for i in range(0,12):
        ax=ke_psi[i,:,:].plot.contourf(x='lon',y='lat',levels=np.arange(0.,250.,10.),extend='both', add_colorbar=False, add_label = False)
        cb1=plt.colorbar(ax)
        cb1.set_label('Nondivergent KE')
        if land_mask:
            land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')
            land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
            land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='0.75',add_colorbar=False,add_labels=False)
        cs = ke_chi[i,:,:].plot.contour(x='lon',y='lat',levels=np.arange(0.,36.,6.), extend='both', add_colorbar=False, colors='k', add_label = False)    
        plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.ylim(-10,45)
        plt.xlim(25,150)
        plt.tight_layout()  
        plt.savefig(plot_dir+str(i+1)+'_'+str(mn_dic[i+1])+'.png')
        plt.close()
        


def abs_vort(run, months, filename='plev_pentad', timeav='month', period_fac=1.,land_mask=False,level=9):
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/abs_vort/'+str(level)+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    mn_dic = month_dic(1)
    
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    abs_vort = data.vor + f
    abs_vort *= 1.e4

    for i in range(0,12):
        ax=abs_vort[i,level,:,:].plot.contourf(x='lon',y='lat',levels=np.arange(-1.,1.1,0.1), extend='both', add_colorbar=False, add_label = False)
        cb1=plt.colorbar(ax)
        cb1.set_label('Absolute vorticity * 10^4, s^{-1}')
        if land_mask:
            land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')
            land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
            land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='0.75',add_colorbar=False,add_labels=False)
        cs = abs_vort[i,level,:,:].plot.contour(x='lon',y='lat',levels=np.arange(-100.,100.,100.), extend='both', add_colorbar=False, colors='k', add_label = False)    
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.ylim(-10,45)
        plt.xlim(25,150)
        plt.tight_layout()  
        plt.savefig(plot_dir+str(i+1)+'_'+str(mn_dic[i+1])+'.png')
        plt.close()


def cross_prod(run, months, filename='plev_pentad', timeav='month', period_fac=1.,land_mask=False,level=9):
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/cross_prod/'+str(level)+'/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    mn_dic = month_dic(1)
    
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    uwnd = data.ucomp[:,level,:,:]
    vwnd = data.vcomp[:,level,:,:]

    w = VectorWind(uwnd, vwnd)
    uchi, vchi, upsi, vpsi = w.helmholtz()
    cross_prod = (upsi*vchi - vpsi*uchi)
    

    for i in range(0,12):
        ax=cross_prod[i,:,:].plot.contourf(x='lon',y='lat',levels=np.arange(-200.,201.,10.),extend='both', add_colorbar=False, add_label = False)
        cb1=plt.colorbar(ax)
        cb1.set_label('Vpsi x Vchi')
        if land_mask:
            land = xr.open_dataset('/scratch/rg419/GFDL_model/GFDLmoistModel/input/land.nc')
            land_plot = xr.DataArray(land.land_mask.values, [('lat', data.lat), ('lon', data.lon)])
            land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='0.75',add_colorbar=False,add_labels=False)
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.ylim(-10,45)
        plt.xlim(25,150)
        plt.tight_layout()  
        plt.savefig(plot_dir+str(i+1)+'_'+str(mn_dic[i+1])+'.png')
        plt.close()        
        
    
sf_vp('ap_1_rd', [121,481])
sf_vp('ap_20', [121,481])
sf_vp('era_amipsst', [121,337],land_mask=True)

ke_components('ap_1_rd', [121,481])
ke_components('ap_20', [121,481])
ke_components('era_amipsst', [121,337],land_mask=True)

abs_vort('ap_1_rd', [121,481])
abs_vort('ap_20', [121,481])
abs_vort('era_amipsst', [121,337],land_mask=True)

cross_prod('ap_1_rd', [121,481])
cross_prod('ap_20', [121,481])
cross_prod('era_amipsst', [121,337],land_mask=True)