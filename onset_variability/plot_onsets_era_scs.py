'''28/09/2018 ERA version of other file: onset dates from 1979 onwards and corresponding Early/Late onset climates in preceeding winter/spring '''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm
from data_handling_updates import gradients as gr
from windspharm.xarray import VectorWind

plot_dir = '/scratch/rg419/plots/onset_variability/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

# Load in onset dates
onsets_scsm = np.load('/scratch/rg419/python_scripts/onset_variability/era_onsets_scsm.npy')
onsets_scsm = xr.DataArray(onsets_scsm, coords={'year': ('year', range(1979,2017))}, dims=['year'])

years = range(1979,2017)

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

# Plot onsets 
plt.plot(years,onsets_scsm,'o-', color='C0')
plt.plot([1979,2016],[round(scsm_stats[0])]*2,'-', color='C0')
plt.fill_between([1979,2016], [round(scsm_stats[0]) - 0.75*scsm_stats[1], round(scsm_stats[0]) - 0.75*scsm_stats[1]], 
                              [round(scsm_stats[0]) + 0.75*scsm_stats[1], round(scsm_stats[0]) + 0.75*scsm_stats[1]], alpha=0.2, color='C0')
plt.xlabel('Year')
plt.ylabel('Onset pentad')
plt.yticks(range(24,37))
plt.grid(True,linestyle=':')
plt.title('South China Sea Monsoon Onset')
plt.savefig(plot_dir + 'scs_onsets_era.pdf', format='pdf')
plt.close()


# Load in ERA data, trim to only use 1979-2016, and resample in months
def load_era_var(var):
    name_temp = '/disca/share/reanalysis_links/monthly_era/' + var + '/interim_' + var + '_%04d.nc'
    names = [name_temp % m for m in years  ]
    data = xr.open_mfdataset(names)
    data = data[var].load().loc['1979-01':'2016-12']
    return data

print('Loading data...')
data_u = load_era_var('u')
data_v = load_era_var('v')
data_vo = load_era_var('vo')
data_q = load_era_var('q')
data_t = load_era_var('t')
data_z = load_era_var('z')

lh = 2.5e6
cp=1005.
data_mse = (data_q * lh + data_t * cp + data_z)/1000.
data_vo = data_vo * 86400.

# Create a VectorWind instance to handle the computation
w = VectorWind(data_u.sel(level=np.arange(50.,950.,50.)), data_v.sel(level=np.arange(50.,950.,50.)))
# Compute variables
streamfun, vel_pot = w.sfvp()
uchi, vchi, upsi, vpsi = w.helmholtz()

coslat = np.cos(data_u.latitude * np.pi/180)

dp=5000.
# Evaluate mass fluxes for the zonal and meridional components (Walker and Hadley) following Schwendike et al. 2014
mass_flux_zon = (gr.ddx(uchi, latname='latitude', lonname='longitude')).cumsum('level') * dp * coslat/ 9.81
mass_flux_merid = (gr.ddy(vchi, latname='latitude')).cumsum('level') * dp * coslat/ 9.81
    
print('Data loaded ok')

land_mask = '/scratch/rg419/python_scripts/land_era/ERA-I_Invariant_0125.nc'
land = xr.open_dataset(land_mask)


def plot_winter_climate(data, title, levels_clim=np.arange(-50.,51.,5.), levels=np.arange(-5.,5.1,0.5), local=False, lev=200.):
    
    # Add a coordinate with the early/late/normal timing    
    data.coords['timing'] = (('time'), np.repeat(onsets_scsm.timing.values,12)) 
    data = data.sel(level=lev)
    
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
    
    f1 = data_clim.plot.contourf(ax=ax1, x='longitude',y='latitude', add_labels=False, add_colorbar=False, extend='both', levels=levels_clim)
    f2 = data_early.plot.contourf(ax=ax2, x='longitude',y='latitude', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    f3 = data_normal.plot.contourf(ax=ax3, x='longitude',y='latitude', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    f4 = data_late.plot.contourf(ax=ax4, x='longitude',y='latitude', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    
    i=0
    for ax in axes:
        land.lsm[0,:,:].plot.contour(ax=ax, x='longitude', y='latitude', levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
        if local:
            ax.set_xlim(70.,150.)
            ax.set_ylim(0.,60.)
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
    
    cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.07, pad=0.2, aspect=30, shrink=1.)
    cb2=fig.colorbar(f2, ax=axes[1:], use_gridspec=True, orientation = 'horizontal',fraction=0.07, pad=0.2, aspect=60, shrink=0.75)
    #cb1.set_label(var)
    
    if local:
        plt.savefig('/scratch/rg419/plots/onset_variability_new/early_vs_late/JFM_' + title + '_era_local.pdf', format='pdf')
    else:
        plt.savefig('/scratch/rg419/plots/onset_variability_new/early_vs_late/JFM_' + title + '_era.pdf', format='pdf')
    plt.close()
    

plot_winter_climate(mass_flux_zon, 'mass_flux_zon_era', levels_clim=np.arange(-0.0055,0.0056,0.001), levels=np.arange(-0.0012,0.0013,0.0002))
plot_winter_climate(mass_flux_zon, 'mass_flux_zon_era', local=True, levels_clim=np.arange(-0.0055,0.0056,0.001), levels=np.arange(-0.0012,0.0013,0.0002))
plot_winter_climate(mass_flux_merid, 'mass_flux_merid_era', levels_clim=np.arange(-0.0055,0.0056,0.001), levels=np.arange(-0.0012,0.0013,0.0002))
plot_winter_climate(mass_flux_merid, 'mass_flux_merid_era', local=True, levels_clim=np.arange(-0.0055,0.0056,0.001), levels=np.arange(-0.0012,0.0013,0.0002))

#plot_winter_climate(data_vo, 'vo', levels_clim=np.arange(-3.,4.,1.), levels=np.arange(-0.6,0.7,0.2))
#plot_winter_climate(data_vo, 'vo', levels_clim=np.arange(-3.,4.,1.), levels=np.arange(-0.6,0.7,0.2), local=True)
#plot_winter_climate(data_u, 'u')
#plot_winter_climate(data_u, 'u', local=True)
#plot_winter_climate(data_mse, 'mse', lev=850., levels_clim=np.arange(250.,340.,10.), levels=np.arange(-2.5,2.6,0.25))
#plot_winter_climate(data_mse, 'mse', local=True, lev=850., levels_clim=np.arange(250.,340.,10.), levels=np.arange(-2.5,2.6,0.25))
#plot_winter_climate(data_t, 't', lev=850., levels_clim=np.arange(240.,301.,10.), levels=np.arange(-1.,1.1,0.1))
#plot_winter_climate(data_t, 't', local=True, lev=850., levels_clim=np.arange(240.,301.,10.), levels=np.arange(-1.,1.1,0.1))
#plot_winter_climate(data_q, 'q', lev=850., levels_clim=np.arange(0,0.015,0.001), levels=np.arange(-0.001,0.0011,0.0002))
#plot_winter_climate(data_q, 'q', local=True, lev=850., levels_clim=np.arange(0,0.015,0.001), levels=np.arange(-0.001,0.0011,0.0002))
#plot_winter_climate(data_z, 'z', lev=850., levels_clim=np.arange(12500.,15001.,200.), levels=np.arange(-100.,101.,10.))
#plot_winter_climate(data_z, 'z', local=True, lev=850., levels_clim=np.arange(12500.,15001.,200.), levels=np.arange(-100.,101.,10.))

data_u.close()
data_vo.close()
data_t.close()
data_q.close()
data_z.close()



