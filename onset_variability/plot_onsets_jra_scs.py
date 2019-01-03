'''26/09/2018 Some figures for Geoff's CSSP talk: onset dates from 1958 onwards and corresponding Early/Late onset climates in preceeding winter/spring '''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm
from data_handling_updates import gradients as gr, model_constants as mc
from windspharm.xarray import VectorWind

plot_dir = '/scratch/rg419/plots/onset_variability/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

# Load in onset dates
onsets_scsm = np.load('/scratch/rg419/python_scripts/onset_variability/jra_onsets_scsm.npy')
onsets_scsm = xr.DataArray(onsets_scsm, coords={'year': ('year', range(1958,2017))}, dims=['year'])

years = range(1958,2017)

# Take mean and standard deviation
scsm_stats = np.mean(np.array(onsets_scsm.sel(year=years))), np.std(np.array(onsets_scsm.sel(year=years)))

# Decide whether monsoon is early, normal, or late
def get_timing(onsets, years):
    onsets = onsets.sel(year=years)
    onset_stats = np.mean(np.array(onsets)), np.std(np.array(onsets))
    enl = []
    enl_no = []
    i=0
    for year in years:
        if onsets[i] > round(onset_stats[0]) + 0.75*onset_stats[1]:
            enl.append('Late')
            enl_no.append(2)
        elif onsets[i] < round(onset_stats[0]) - 0.75*onset_stats[1]:
            enl.append('Early')
            enl_no.append(0)
        else:
            enl.append('Normal')
            enl_no.append(1)
            
        i=i+1    
    onsets_out = xr.DataArray(onsets, coords={'timing': ('year', enl), 'timing_no': ('year', enl_no), 'year': ('year', years)}, dims=['year'])
    return onsets_out
#print(onsets_scsm.swap_dims({'year': 'timing'}).sel(timing='Normal').year)

onsets_scsm = get_timing(onsets_scsm, years)
#print(onsets_scsm.timing.values)

# Plot onsets 
plt.plot(years,onsets_scsm,'o-', color='C0')
plt.plot([1958,2016],[round(scsm_stats[0])]*2,'-', color='C0')
plt.fill_between([1958,2016], [round(scsm_stats[0]) - 0.75*scsm_stats[1], round(scsm_stats[0]) - 0.75*scsm_stats[1]], 
                              [round(scsm_stats[0]) + 0.75*scsm_stats[1], round(scsm_stats[0]) + 0.75*scsm_stats[1]], alpha=0.2, color='C0')
plt.xlabel('Year')
plt.ylabel('Onset pentad')
plt.yticks(range(24,38,2))
plt.grid(True,linestyle=':')
plt.title('South China Sea Monsoon Onset')
plt.savefig(plot_dir + 'scs_onsets_jra.pdf', format='pdf')
plt.close()



# Load in JRA data
data_u = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/ucomp_monthly/atmos_monthly_together.nc')
data_u = data_u['var33'].load().loc['1958-01':'2016-12']

# v has different time coord to u, presumably due to how Stephen has downloaded/averaged. I think the two are equivalent, so just substitute the time dimension into v
data_v_temp = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/vcomp_monthly/atmos_monthly_together.nc')
data_v = xr.DataArray(data_v_temp.var34.values, coords={'time': data_u.time, 'lev': data_v_temp.lev,
                                                  'lat': data_v_temp.lat, 'lon': data_v_temp.lon}, dims=['time','lev','lat','lon'])

dvdx = gr.ddx(data_v, a=6371.0e3) 
dudy = gr.ddy(data_u, a=6371.0e3)
data_vo = (dvdx - dudy)*86400.

data_t = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/temp_monthly/atmos_monthly_together.nc')
data_q = xr.open_dataset('/disca/share/rg419/JRA_55/sphum_monthly/atmos_monthly_together.nc')
data_z = xr.open_dataset('/disca/share/reanalysis_links/jra_55/1958_2016/height_monthly/atmos_monthly_together.nc')

data_mse = (mc.cp_air*data_t.var11 + mc.L*data_q.var51 + 9.81*data_z.var7)/1000.

# Create a VectorWind instance to handle the computation
w = VectorWind(data_u.sel(lev=np.arange(5000.,100001.,5000.)), data_v.sel(lev=np.arange(5000.,100001.,5000.)))
# Compute variables
streamfun, vel_pot = w.sfvp()
uchi, vchi, upsi, vpsi = w.helmholtz()
coslat = np.cos(data_u.lat * np.pi/180)

dp=5000.
# Evaluate mass fluxes for the zonal and meridional components (Walker and Hadley) following Schwendike et al. 2014
mass_flux_zon = (gr.ddx(uchi)).cumsum('lev') * dp * coslat/ 9.81
mass_flux_merid = (gr.ddy(vchi)).cumsum('lev') * dp * coslat/ 9.81

land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
land = xr.open_dataset(land_mask)


def plot_winter_climate(data, title, levels_clim=np.arange(-50.,51.,5.), levels=np.arange(-5.,5.1,0.5), local=False, lev=20000.):
    
    # Add a coordinate with the early/late/normal timing    
    data.coords['timing'] = (('time'), np.repeat(onsets_scsm.timing.values,12)) 
    data = data.sel(lev=lev)
    
    data_clim = data.groupby('time.month').mean('time').sel(month=[1,2,3]).mean('month')
    data_early = (data.where(data['timing']=='Early', drop=True).groupby('time.month').mean('time') - data.groupby('time.month').mean('time')).sel(month=[1,2,3]).mean('month')
    data_normal = (data.where(data['timing']=='Normal', drop=True).groupby('time.month').mean('time') - data.groupby('time.month').mean('time')).sel(month=[1,2,3]).mean('month')
    data_late = (data.where(data['timing']=='Late', drop=True).groupby('time.month').mean('time') - data.groupby('time.month').mean('time')).sel(month=[1,2,3]).mean('month')
    
    # Start figure with 4 subplots
    rcParams['figure.figsize'] = 15, 4.5
    rcParams['font.size'] = 14
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey='row')
    axes = [ax1, ax2, ax3, ax4]
    title_list = ['Climatology', 'Early onset', 'Normal onset', 'Late onset']
    
    f1 = data_clim.plot.contourf(ax=ax1, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels_clim)
    f2 = data_early.plot.contourf(ax=ax2, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    f3 = data_normal.plot.contourf(ax=ax3, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    f4 = data_late.plot.contourf(ax=ax4, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    
    i=0
    for ax in axes:
        land.lsm[0,:,:].plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        if local:
            ax.set_xlim(70.,150.)
            ax.set_ylim(-10.,50.)
        else:
            ax.set_xticks(np.arange(0.,361.,90.))
            ax.set_yticks(np.arange(-60.,61.,30.))
            ax.set_ylim(-60.,60.)
        ax.set_title(title_list[i])
        ax.grid(True,linestyle=':')
        ax.set_xlabel('Longitude')
        i=i+1
    
    ax1.set_ylabel('Latitude')
    
    plt.subplots_adjust(left=0.06, right=0.97, top=0.92, bottom=0.1, hspace=0.3, wspace=0.2)
    
    
    if title=='mse_jra':
        cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.07, pad=0.2, aspect=60, shrink=1., ticks=np.arange(250.,340.,20.))
    else:
        cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.07, pad=0.2, aspect=30, shrink=1.)
    cb2=fig.colorbar(f2, ax=axes[1:], use_gridspec=True, orientation = 'horizontal',fraction=0.07, pad=0.2, aspect=60, shrink=0.75)
    #cb1.set_label(var)
    
    if local:
        plt.savefig('/scratch/rg419/plots/onset_variability_new/early_vs_late/JFM_' + title + '_local.pdf', format='pdf')
    else:
        plt.savefig('/scratch/rg419/plots/onset_variability_new/early_vs_late/JFM_' + title + '.pdf', format='pdf')
    plt.close()
    

plot_winter_climate(data_vo, 'vo_jra', levels_clim=np.arange(-3.,4.,1.), levels=np.arange(-0.6,0.7,0.2))
plot_winter_climate(data_vo, 'vo_jra', levels_clim=np.arange(-3.,4.,1.), levels=np.arange(-0.6,0.7,0.2), local=True)
plot_winter_climate(data_u, 'ucomp')
plot_winter_climate(data_u, 'ucomp', local=True)
plot_winter_climate(data_mse, 'mse', local=True, levels_clim=np.arange(250.,331.,10.), levels=np.arange(-2.5,2.7,0.25), lev=85000.)
plot_winter_climate(data_mse, 'mse', levels_clim=np.arange(250.,331.,10.), levels=np.arange(-2.5,2.7,0.25), lev=85000.)
#plot_winter_climate(mass_flux_zon, 'mass_flux_zon_jra', levels_clim=np.arange(-0.0055,0.0056,0.001), levels=np.arange(-0.0012,0.0013,0.0002))
#plot_winter_climate(mass_flux_zon, 'mass_flux_zon_jra', local=True, levels_clim=np.arange(-0.0055,0.0056,0.001), levels=np.arange(-0.0012,0.0013,0.0002))
#plot_winter_climate(mass_flux_merid, 'mass_flux_merid_jra', levels_clim=np.arange(-0.0055,0.0056,0.001), levels=np.arange(-0.0012,0.0013,0.0002))
#plot_winter_climate(mass_flux_merid, 'mass_flux_merid_jra', local=True, levels_clim=np.arange(-0.0055,0.0056,0.001), levels=np.arange(-0.0012,0.0013,0.0002))

data_u.close()
data_v.close()
data_vo.close()



