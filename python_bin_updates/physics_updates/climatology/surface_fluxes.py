# 1/12/2017 Plot hms of SST, precip, SW net, LW net, SH and LH

from data_handling_updates import month_dic,  model_constants as mc, gradients as gr, make_sym
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from climatology import precip_centroid
from hadley_cell import mass_streamfunction

g = 9.8
cp = 287.04/2*7
L = 2.500e6
rho_cp = 1.035e3 * 3989.24495292815

def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons


def surface_plot(run, lonin=[-1.,361.], diff_run=None, scale_fac=1., do_make_sym=True):
    
    print('what?')
    rcParams['figure.figsize'] = 15, 10
    rcParams['font.size'] = 14
    
    plot_dir = '/scratch/rg419/plots/surface_fluxes/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    #Calculate psi to overplot
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    psi = make_sym(psi, asym=True)
    
    # Make precip symmetric and find the precip centroid
    data['precipitation'] = make_sym(data.precipitation)
    precip_centroid(data)
    
    data['flux_lw_up'] = data.t_surf ** 4. * mc.stefan
    data['dTdt'] = gr.ddt(data.t_surf) * 86400. * scale_fac
    
    if do_make_sym:
        for field in ['t_surf','dTdt','flux_sw','flux_lw','flux_lw_up','flux_t','flux_lhe']:
            data[field] = make_sym(data[field])
    
    lons = pick_lons(data, lonin)
    
    if not diff_run == None:
        diff_data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + diff_run + '.nc')
        data = (data - diff_data).sel(lon=lons).mean('lon')
        levels_t = np.arange(-10.,11.,1.)
        levels_dt = np.arange(-0.5,0.51,0.05)
        levels_flux = np.arange(-80.,81.,5.)
    else:
        data = data.sel(lon=lons).mean('lon')
        levels_t = np.arange(270.,306.,2.5)
        levels_dt = np.arange(-0.3,0.31,0.05)
        levels_flux = np.arange(-300.,305.,20.)
    

    
    f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
    left_column = [ax1,ax3,ax5]
    bottom_row = [ax5,ax6]
    all_plots = [ax1,ax2,ax3,ax4,ax5,ax6]
    
    #check = data.flux_sw + (data.flux_lw - flux_lw_up) - data.flux_t - data.flux_lhe
    #check = data.flux_sw - data.flux_lhe
    
    data.t_surf.plot.contourf(x='xofyear', y='lat', levels=levels_t, ax=ax1, extend = 'both', add_labels=False)
    ax1.set_title('SST')
    
    data.dTdt.plot.contourf(x='xofyear', y='lat',  ax=ax2, levels=levels_dt, extend = 'both', add_labels=False, cmap='RdBu_r')
    ax2.set_title('d(SST)/dt')
    
    data.flux_sw.plot.contourf(x='xofyear', y='lat', levels=levels_flux, ax=ax3, extend = 'both', add_labels=False, cmap='RdBu_r')
    ax3.set_title('Net downward SW flux')
    
    (data.flux_lw - data.flux_lw_up).plot.contourf(x='xofyear', y='lat', levels=levels_flux, ax=ax4, extend = 'both', add_labels=False, cmap='RdBu_r')
    ax4.set_title('Net downward LW flux')
    
    (-1.*data.flux_t).plot.contourf(x='xofyear', y='lat', levels=levels_flux, ax=ax5, extend = 'both', add_labels=False, cmap='RdBu_r')
    #check.plot.contourf(x='xofyear', y='lat', levels=np.arange(-300.,305.,20.), ax=ax5, extend = 'both', add_labels=False, cmap='RdBu_r')
    ax5.set_title('Downward sensible heat flux')
    
    (-1.*data.flux_lhe).plot.contourf(x='xofyear', y='lat', levels=levels_flux, ax=ax6, extend = 'both', add_labels=False, cmap='RdBu_r')
    ax6.set_title('Downward latent heat flux')
    
    
    i=0
    labels=['a)','b)','c)','d)','e)','f)']
    for ax in all_plots:
        ax.grid(True,linestyle=':')
        ax.set_xticks([12,24,36,48,60,72])
        ax.set_ylim([-60,60])
        ax.set_yticks([-60,-30,0,30,60])
        ax.text(-8, 60., labels[i])
        psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
        psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
        psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
        data.p_cent.plot.line(ax=ax, color='k', linewidth=2)
        ax.set_ylabel('')
        ax.set_xlabel('')
        i=i+1
    
    for ax in left_column:
        ax.set_ylabel('Latitude')
    for ax in bottom_row:
        ax.set_xlabel('Pentad')
        
    plt.subplots_adjust(left=0.07, right=0.97, top=0.97, bottom=0.05, hspace=0.2, wspace=0.1)
    
    if lonin == [-1.,361.]:
        if not diff_run == None:
            plt.savefig(plot_dir + 'surface_fluxes_' + run + '_' + diff_run + '.pdf', format='pdf')
        
        else:
            plt.savefig(plot_dir + 'surface_fluxes_' + run + '.pdf', format='pdf')
    else:
        figname = plot_dir + 'surface_fluxes_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'
        plt.savefig(figname, format='pdf')
        
    plt.close()



def surf_cooling_plot(run, lonin=[-1.,361.], mld=10.):
    
    rcParams['figure.figsize'] = 6, 10
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/surface_fluxes/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    lons = pick_lons(data, lonin)
    
    data = data.sel(lon=lons).mean('lon')
    
    flux_lw_up = data.t_surf ** 4. * mc.stefan
    dTdt = gr.ddt(data.t_surf) * 86400.
    
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    
    sum_of_all = (data.flux_sw + (data.flux_lw - flux_lw_up) - data.flux_t - data.flux_lhe)/rho_cp/mld* 86400. 
    sw_and_lhe = (data.flux_sw - data.flux_lhe + (data.flux_lw - flux_lw_up - data.flux_t).mean('xofyear'))/rho_cp/mld* 86400. 
    sw_and_lw = (data.flux_sw + (data.flux_lw - flux_lw_up) - (data.flux_lhe + data.flux_t).mean('xofyear'))/rho_cp/mld* 86400. 
    
    f1 = sum_of_all.plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.,1.05,0.1), ax=ax1, extend = 'both', add_labels=False, add_colorbar=False)
    ax1.set_title('Sum of all fluxes')
    
    sw_and_lhe.plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.,1.05,0.1), ax=ax2,  extend = 'both', add_labels=False, cmap='RdBu_r', add_colorbar=False)
    ax2.set_title('SW + LHE')
    
    sw_and_lw.plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.,1.05,0.1), ax=ax3, extend = 'both', add_labels=False, cmap='RdBu_r', add_colorbar=False)
    ax3.set_title('SW + LW')
    
    ax3.set_xlabel('Pentad')
    for ax in [ax1,ax2,ax3]:
        ax.set_ylabel('Latitude')
        ax.grid(True,linestyle=':')
    
    plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
    #Colorbar
    cb1=f.colorbar(f1, ax=(ax1, ax2, ax3), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
    cb1.set_label('Downward flux, W/m^2')
    
    plt.savefig(plot_dir + 'surf_cooling_plot_' + run + '.pdf', format='pdf')
    plt.close()



def int_flux_plot(run, lonin=[-1.,361.], mld=10.):
    
    rcParams['figure.figsize'] = 6, 13.3
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/surface_fluxes/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    lons = pick_lons(data, lonin)
    
    data = data.sel(lon=lons).mean('lon')
    
    flux_lw_up = data.t_surf ** 4. * mc.stefan
    dTdt = gr.ddt(data.t_surf) * 86400.
    
    
    def int_flux_fwd(flux):
        #cp dTdt = SW + LW + LH + SH
        n = data.t_surf.shape[0]
        t_out = np.zeros([n,64])
        t_out[0,:] = data.t_surf[0,:]   
        for i in range(n-1):
            t_out[i+1,:] = t_out[i,:] + flux[i,:]/rho_cp/mld * 86400. * 5.
        t_out = xr.DataArray(t_out, coords=[data.xofyear.values, data.lat.values], dims=['xofyear', 'lat'])
        return t_out
    
    t_all = int_flux_fwd(data.flux_sw - data.flux_lhe + data.flux_lw - flux_lw_up - data.flux_t)
    t_sw = int_flux_fwd(data.flux_sw + (-data.flux_lhe + data.flux_lw - flux_lw_up - data.flux_t).mean('xofyear'))
    t_lw = int_flux_fwd(data.flux_sw + data.flux_lw - flux_lw_up - data.flux_lhe.mean('xofyear'))
    t_lh = int_flux_fwd(data.flux_sw - data.flux_lhe + (data.flux_lw - flux_lw_up).mean('xofyear'))
    
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)
    
    
    f1 = t_all.plot.contourf(x='xofyear', y='lat', levels=np.arange(250.,311.,5.), ax=ax1, extend = 'both', add_labels=False, add_colorbar=False)
    ax1.set_title('Sum of all fluxes')
    
    f1 = t_sw.plot.contourf(x='xofyear', y='lat', levels=np.arange(250.,311.,5.), ax=ax2, extend = 'both', add_labels=False, add_colorbar=False)
    ax2.set_title('SW')
    
    t_lw.plot.contourf(x='xofyear', y='lat', levels=np.arange(250.,311.,5.),  ax=ax3,  extend = 'both', add_labels=False, add_colorbar=False)
    ax3.set_title('SW + LHE')
    
    t_lh.plot.contourf(x='xofyear', y='lat', levels=np.arange(250.,311.,5.),  ax=ax4, extend = 'both', add_labels=False, add_colorbar=False)
    ax4.set_title('SW + LW')
    
    ax4.set_xlabel('Pentad')
    for ax in [ax1,ax2,ax3, ax4]:
        ax.set_ylabel('Latitude')
        ax.grid(True,linestyle=':')
    
    plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
    #Colorbar
    cb1=f.colorbar(f1, ax=(ax1, ax2, ax3, ax4), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
    cb1.set_label('Temperature, K')
    
    plt.savefig(plot_dir + 'int_flux_plot_' + run + '.pdf', format='pdf')
    plt.close()
    


def fixed_sst_imbalance(run, lonin=[-1.,361.], mld=10.):
    # Plot the imbalance between actual rate of change of temperature and that the fluxes should cause for fixed SST expts
    
    rcParams['figure.figsize'] = 6, 10
    rcParams['font.size'] = 16
    
    plot_dir = '/scratch/rg419/plots/surface_fluxes/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    lons = pick_lons(data, lonin)
    
    data = data.sel(lon=lons).mean('lon')
    
    flux_lw_up = data.t_surf ** 4. * mc.stefan
    dTdt = gr.ddt(data.t_surf) * 86400.
    
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    
    sum_of_all = (data.flux_sw + (data.flux_lw - flux_lw_up) - data.flux_t - data.flux_lhe)/rho_cp/mld* 86400. 
    dTdt = gr.ddt(data.t_surf) * 86400.

    
    f1 = sum_of_all.plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.,1.,0.1), ax=ax1, extend = 'both', add_labels=False, add_colorbar=False)
    ax1.set_title('Sum of all fluxes')
    
    dTdt.plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.,1.,0.1), ax=ax2,  extend = 'both', add_labels=False, cmap='RdBu_r', add_colorbar=False)
    ax2.set_title('dTs/dt')
    
    (sum_of_all - dTdt).plot.contourf(x='xofyear', y='lat', levels=np.arange(-1.,1.,0.1), ax=ax3, extend = 'both', add_labels=False, cmap='RdBu_r', add_colorbar=False)
    ax3.set_title('Fluxes-dTsdt')
    
    ax3.set_xlabel('Pentad')
    for ax in [ax1,ax2,ax3]:
        ax.set_ylabel('Latitude')
        ax.grid(True,linestyle=':')
    
    plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
    #Colorbar
    cb1=f.colorbar(f1, ax=(ax1, ax2, ax3), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
    cb1.set_label('dT/dt, K/day')
    
    plt.savefig(plot_dir + 'fixed_sst_imbalance_' + run + '.pdf', format='pdf')
    plt.close()
    
    
    

if __name__ == "__main__":
    
    #for run in ['sn_1.000', 'sn_0.500', 'sn_2.000', 'sine_sst_10m', 'sine_sst_10m_zs', 'rt_0.500', 'rt_2.000', 'zs_sst']:
    #    print run
    #    surf_cooling_plot(run)
    #surf_cooling_plot('ap_2', mld=2.)
    #surf_cooling_plot('ap_20', mld=20.)
    surface_plot('half_shallow', lonin=[340,20])
    surface_plot('half_shallow', lonin=[70,110])
    surface_plot('half_shallow', lonin=[160,200])
    surface_plot('half_shallow', lonin=[250,290])

    
