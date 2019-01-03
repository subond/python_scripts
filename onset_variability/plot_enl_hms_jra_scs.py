'''04/12/2018 Plot averaged hms for early, normal and lat SCS monsoon onset'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm
from data_handling_updates import gradients as gr, model_constants as mc
from windspharm.xarray import VectorWind


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


def pentad_means_of_year(data, years):  # Function to get pentad of year
    data_pentad = np.zeros( (len(years), 73, len(data.lat), len(data.lon)) )
    i=0
    for year in years:
        data_year = data.sel(time=str(year))
        if len(data_year.time)==366:
            pentad = np.repeat(np.arange(1., 74.), 5)
            pentad = np.insert(pentad, 10, 2)    
        else:
            pentad = np.repeat(np.arange(1., 74.), 5)
        data_year = data_year.assign_coords(pentad = ('time', pentad))
        data_year = data_year.groupby('pentad').mean(('time'))
        data_pentad[i,:,:,:] = data_year.values
        i=i+1
    data_pentad = xr.DataArray(data_pentad, 
           coords={'year': ('year', years), 'pentad': ('pentad', data_year.pentad), 'lat': ('lat', data.lat), 'lon':('lon', data.lon)}, 
           dims=['year', 'pentad', 'lat', 'lon'])
    return data_pentad
        
    

def timing_hms(data, title, lonin=[110.,120.], levels=np.arange(-50.,51.,5.), anom=True):
    
    data = pentad_means_of_year(data, np.arange(1958,2017))
    data.coords['timing'] = (('year'), onsets_scsm.timing)
    
    lons = [data_u.lon[i].values for i in range(len(data_u.lon)) if data_u.lon[i] >= lonin[0] and data_u.lon[i] <= lonin[1]]
    
    data = data.sel(lon=lons).mean('lon')
    
    if anom:
        data_early = (data.where(data['timing']=='Early', drop=True).mean('year') - data.mean('year'))
        data_normal = (data.where(data['timing']=='Normal', drop=True).mean('year') - data.mean('year'))
        data_late = (data.where(data['timing']=='Late', drop=True).mean('year') - data.mean('year'))
    else:
        data_early = data.where(data['timing']=='Early', drop=True).mean('year')
        data_normal = data.where(data['timing']=='Normal', drop=True).mean('year')
        data_late = data.where(data['timing']=='Late', drop=True).mean('year')
    
    # Start figure with 3 subplots
    rcParams['figure.figsize'] = 15, 4.5
    rcParams['font.size'] = 14
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey='row')
    axes = [ax1, ax2, ax3]
    title_list = ['Early onset', 'Normal onset', 'Late onset']
    
    f1 = data_early.plot.contourf(ax=ax1, x='pentad',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    f2 = data_normal.plot.contourf(ax=ax2, x='pentad',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    f3 = data_late.plot.contourf(ax=ax3, x='pentad',y='lat', add_labels=False, add_colorbar=False, extend='both', levels=levels)
    
    i=0
    for ax in axes:
        ax.set_title(title_list[i])
        ax.grid(True,linestyle=':')
        ax.set_xlim(0.,73.)
        ax.set_xlabel('Pentad')
        i=i+1
    
    ax1.set_ylabel('Latitude')
    
    plt.subplots_adjust(left=0.06, right=0.97, top=0.92, bottom=0.1, hspace=0.3, wspace=0.2)
    
    cb2=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.07, pad=0.2, aspect=60, shrink=0.75)
    #cb1.set_label(var)
    
    if anom:
        plt.savefig('/scratch/rg419/plots/onset_variability_new/early_vs_late/hm_' + title + '_anom.pdf', format='pdf')
    else:
        plt.savefig('/scratch/rg419/plots/onset_variability_new/early_vs_late/hm_' + title + '.pdf', format='pdf')
    plt.close()


plot_dir = '/scratch/rg419/plots/onset_variability_new/early_vs_late/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
years = range(1958,2017)

# Load in onset dates
onsets_scsm = np.load('/scratch/rg419/python_scripts/onset_variability/jra_onsets_scsm.npy')
onsets_scsm = xr.DataArray(onsets_scsm, coords={'year': ('year', years)}, dims=['year'])
onsets_scsm = get_timing(onsets_scsm, years)

# Load in JRA data
data_u = xr.open_dataset('/disca/share/rg419/jra_ucomp_daily_200.nc', chunks={'time': 30})
data_u = data_u['var33'].load().loc['1958-01':'2016-12']
# v has different time coord to u, presumably due to how Stephen has downloaded/averaged. I think the two are equivalent, so just substitute the time dimension into v
data_v_temp = xr.open_dataset('/disca/share/rg419/jra_vcomp_daily_200.nc', chunks={'time': 30})
data_v = xr.DataArray(data_v_temp.var34.values, coords={'time': data_u.time, 'lat': data_u.lat, 'lon': data_u.lon}, dims=('time','lat','lon'))    
print('files opened')

dvdx = gr.ddx(data_v, a=6371.0e3) 
dudy = gr.ddy(data_u, a=6371.0e3)
data_vo = (dvdx - dudy)*86400.

timing_hms(data_vo, 'vo_jra', levels=np.arange(-3.,4.,1.), anom=False)
timing_hms(data_u, 'ucomp', anom=False)
timing_hms(data_vo, 'vo_jra', levels=np.arange(-0.6,0.7,0.2))
timing_hms(data_u, 'ucomp', levels=np.arange(-5.,5.1,0.5))

data_u.close()
data_v.close()
data_vo.close()

data_t = xr.open_dataset('/disca/share/rg419/jra_temp_daily_850.nc', chunks={'time': 30})
data_q = xr.open_dataset('/disca/share/rg419/jra_sphum_daily_850.nc', chunks={'time': 30})
data_z = xr.open_dataset('/disca/share/rg419/jra_height_daily_850.nc', chunks={'time': 30})

data_mse = (mc.cp_air*data_t.var11 + mc.L*data_q.var51 + 9.81*data_z.var7).sel(lev=85000.)/1000.

timing_hms(data_mse, 'mse_jra', levels=np.arange(250.,331.,10.), anom=False)
timing_hms(data_mse, 'mse_jra', levels=np.arange(-2.5,2.7,0.25))



