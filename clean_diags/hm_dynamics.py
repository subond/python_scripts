"""
Functions to plot Hovmoller at 850hPa of
1) Cross product of irrotational and nondivergent wind
2) Kinetic energy associated with irrotational and nondivergent motions
3) Absolute vorticity

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from data_handling import time_means, cell_area, month_dic
from windspharm.xarray import VectorWind
import sh


def cross_prod(run, months, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.],level=10):
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/hm/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    mn_dic = month_dic(1)
    tickspace = range(7,72,12)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    if lonin[1]>lonin[0]:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]

    
    uwnd = data.ucomp[:,level,:,:].load()
    vwnd = data.vcomp[:,level,:,:].load()

    w = VectorWind(uwnd, vwnd)
    uchi, vchi, upsi, vpsi = w.helmholtz()
    cross_prod = (upsi*vchi - vpsi*uchi)[:,:,lons].mean('lon')
    
    ax=cross_prod.plot.contourf(x='xofyear',y='lat',levels=np.arange(-200.,201.,10.), extend='both', add_colorbar=False, add_label = False)
    uwnd[:,:,lons].mean('lon').plot.contour(x='xofyear',y='lat',levels=np.arange(-1000.,1001.,1000.), colors='k', add_colorbar=False, add_label = False)    
    plt.grid(True,linestyle=':')
    cb1=plt.colorbar(ax)
    cb1.set_label('Vpsi x Vchi')
    plt.ylabel('Latitude')
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.yticks(range(-60,61,30))
    plt.tight_layout()  
    plt.savefig(plot_dir+'cross_prod.png')
    plt.close()



def ke_hm(run, months, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.],level=10):
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/hm/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    mn_dic = month_dic(1)
    tickspace = range(7,72,12)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    if lonin[1]>lonin[0]:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        
    uwnd = data.ucomp[:,level,:,:].load()
    vwnd = data.vcomp[:,level,:,:].load()

    w = VectorWind(uwnd, vwnd)
    uchi, vchi, upsi, vpsi = w.helmholtz()
    ke_chi = 0.5*(uchi*uchi + vchi*vchi)[:,:,lons].mean('lon')
    ke_psi = 0.5*(upsi*upsi + vpsi*vpsi)[:,:,lons].mean('lon')
    
    ax=ke_psi.plot.contourf(x='xofyear',y='lat',levels=np.arange(0.,250.,10.),extend='both', add_colorbar=False, add_label = False)
    plt.grid(True,linestyle=':')
    cb1=plt.colorbar(ax)
    cb1.set_label('Nondivergent KE')
    cs = uwnd[:,:,lons].mean('lon').plot.contour(x='xofyear',y='lat',levels=np.arange(-1000.,1001.,1000.), extend='both', colors='w', add_colorbar=False, add_label = False)    
    cs = ke_chi.plot.contour(x='xofyear',y='lat',levels=np.arange(0.,36.,6.), extend='both', add_colorbar=False, colors='k', add_label = False)    
    plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.yticks(range(-60,61,30))
    plt.ylabel('Latitude')
    plt.tight_layout()  
    plt.savefig(plot_dir+'ke.png')
    plt.close()
        


def abs_vort_hm(run, months, filename='plev_pentad', timeav='pentad', period_fac=1.,lonin=[-1.,361.],level=10):
    
    plot_dir = '/scratch/rg419/plots/clean_diags/'+run+'/hm/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    mn_dic = month_dic(1)
    tickspace = range(7,72,12)
    labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
    data = time_means(run, months, filename=filename, timeav=timeav,period_fac=period_fac)
    
    if lonin[1]>lonin[0]:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [i for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        
    omega = 7.2921150e-5
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    abs_vort = (data.vor + f)[:,:,:,lons].mean('lon')
    abs_vort *= 1.e4
    
    uwnd = data.ucomp[:,level,:,:].load()

    ax=abs_vort[:,level,:].plot.contourf(x='xofyear',y='lat',levels=np.arange(-1.,1.1,0.1),extend='both', add_colorbar=False, add_label = False)
    cs = uwnd[:,:,lons].mean('lon').plot.contour(x='xofyear',y='lat',levels=np.arange(-1000.,1001.,1000.), extend='both', colors='k', add_colorbar=False, add_label = False)    
    plt.grid(True,linestyle=':')
    cb1=plt.colorbar(ax)
    cb1.set_label('Absolute vorticity')
    plt.xlabel('')
    plt.xticks(tickspace,labels,rotation=25)
    plt.yticks(range(-60,61,30))
    plt.ylabel('Latitude')
    plt.tight_layout()  
    plt.savefig(plot_dir+'abs_vort.png')
    plt.close()
        
        
    
cross_prod('ap_1_rd', [121,481])
cross_prod('ap_20', [121,481])
cross_prod('era_amipsst', [121,337],lonin=[25.,150.])

ke_hm('ap_1_rd', [121,481])
ke_hm('ap_20', [121,481])
ke_hm('era_amipsst', [121,337],lonin=[25.,150.])

abs_vort_hm('ap_1_rd', [121,481])
abs_vort_hm('ap_20', [121,481])
abs_vort_hm('era_amipsst', [121,337],lonin=[25.,150.])

