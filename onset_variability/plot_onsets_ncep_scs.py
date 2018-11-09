'''28/09/2018 Validate onset dates for SCS monsoon - compare NCEP-NCAR with Wang et al. 2004 '''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm
from data_handling_updates import gradients as gr

plot_dir = '/scratch/rg419/plots/onset_variability/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

# Load in onset dates
onsets_scsm = np.load('/scratch/rg419/python_scripts/onset_variability/ncep_onsets_scsm.npy')
onsets_scsm = xr.DataArray(onsets_scsm, coords={'year': ('year', range(1948,2017))}, dims=['year'])

wang_onset = np.array([26,30,26,25,27,26,31,29,31,32,29,30,30,27,28,30,28,29,25,29,34,29,32,25,26,33,29,
                       31,26,28,29,27,27,31,31,31,29,30,27,32,29,28,28,32,28,32,25,27,26,28,29,30,26,26])

wang_onset = xr.DataArray(wang_onset, coords={'year': ('year', range(1948,2002))}, dims=['year'])

years = range(1948,2017)

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
plt.plot(wang_onset.year,wang_onset,'o:', markerfacecolor='none', color='k')
plt.plot([1948,2016],[round(scsm_stats[0])]*2,'-', color='C0')
plt.fill_between([1948,2016], [round(scsm_stats[0]) - 0.75*scsm_stats[1], round(scsm_stats[0]) - 0.75*scsm_stats[1]], 
                              [round(scsm_stats[0]) + 0.75*scsm_stats[1], round(scsm_stats[0]) + 0.75*scsm_stats[1]], alpha=0.2, color='C0')
plt.xlabel('Year')
plt.ylabel('Onset pentad')
plt.yticks(range(24,37))
plt.grid(True,linestyle=':')
plt.title('South China Sea Monsoon Onset')
plt.savefig(plot_dir + 'scs_onsets_ncep.pdf', format='pdf')
plt.close()


# Load in NCEP data, trim to only use 1948-2016
data_u = xr.open_mfdataset('/disca/share/reanalysis_links/NCEP_NCAR/uwnd.mon.mean.nc')
data_u = data_u['uwnd'].load().loc['1948-01':'2016-12']
print('u loaded')
data_v = xr.open_mfdataset('/disca/share/reanalysis_links/NCEP_NCAR/vwnd.mon.mean.nc')
data_v = data_v['vwnd'].load().loc['1948-01':'2016-12']
print('v loaded')
data_t = xr.open_mfdataset('/disca/share/reanalysis_links/NCEP_NCAR/air.mon.mean.nc')
data_t = data_t['air'].load().loc['1948-01':'2016-12'] + 273.15
print('t loaded')
data_q = xr.open_mfdataset('/disca/share/reanalysis_links/NCEP_NCAR/shum.mon.mean.nc')
data_q = data_q['shum'].load().loc['1948-01':'2016-12'] /1000.
print('q loaded')
data_z = xr.open_mfdataset('/disca/share/reanalysis_links/NCEP_NCAR/hgt.mon.mean.nc')
data_z = data_z['hgt'].load().loc['1948-01':'2016-12'] * 9.81
print('z loaded')

dvdx = gr.ddx(data_v, a=6371.0e3) 
dudy = gr.ddy(data_u, a=6371.0e3)
data_vo = (dvdx.values - dudy.values)*86400.
data_vo = xr.DataArray( data_vo, dims = data_u.dims, coords = data_u.coords )
lh = 2.5e6
cp=1005.
data_mse = (data_q * lh + data_t * cp + data_z)/1000.
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
    
    f1 = data_clim.plot.contourf(ax=ax1, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels_clim)
    f2 = data_early.plot.contourf(ax=ax2, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    f3 = data_normal.plot.contourf(ax=ax3, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    f4 = data_late.plot.contourf(ax=ax4, x='lon',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    
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
        plt.savefig('/scratch/rg419/plots/onset_variability_new/early_vs_late/JFM_' + title + '_ncep_local.pdf', format='pdf')
    else:
        plt.savefig('/scratch/rg419/plots/onset_variability_new/early_vs_late/JFM_' + title + '_ncep.pdf', format='pdf')
    plt.close()
    

plot_winter_climate(data_vo, 'vo', levels_clim=np.arange(-3.,4.,1.), levels=np.arange(-0.6,0.7,0.2))
plot_winter_climate(data_vo, 'vo', levels_clim=np.arange(-3.,4.,1.), levels=np.arange(-0.6,0.7,0.2), local=True)
plot_winter_climate(data_u, 'u')
plot_winter_climate(data_u, 'u', local=True)
plot_winter_climate(data_mse, 'mse', lev=850., levels_clim=np.arange(250.,340.,10.), levels=np.arange(-2.5,2.6,0.25))
plot_winter_climate(data_mse, 'mse', local=True, lev=850., levels_clim=np.arange(250.,340.,10.), levels=np.arange(-2.5,2.6,0.25))
plot_winter_climate(data_t, 't', lev=850., levels_clim=np.arange(240.,301.,10.), levels=np.arange(-1.,1.1,0.1))
plot_winter_climate(data_t, 't', local=True, lev=850., levels_clim=np.arange(240.,301.,10.), levels=np.arange(-1.,1.1,0.1))
plot_winter_climate(data_q, 'q', lev=850., levels_clim=np.arange(0,0.015,0.001), levels=np.arange(-0.001,0.0011,0.0002))
plot_winter_climate(data_q, 'q', local=True, lev=850., levels_clim=np.arange(0,0.015,0.001), levels=np.arange(-0.001,0.0011,0.0002))
plot_winter_climate(data_z, 'z', lev=850., levels_clim=np.arange(12500.,15001.,200.), levels=np.arange(-100.,101.,10.))
plot_winter_climate(data_z, 'z', local=True, lev=850., levels_clim=np.arange(12500.,15001.,200.), levels=np.arange(-100.,101.,10.))

data_u.close()
data_vo.close()
data_t.close()
data_q.close()
data_z.close()




