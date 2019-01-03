''' 
7/12/2018 Make Isca hms of different regions in idealised asymmetry runs e.g. half shallow
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from data_handling_updates import month_dic, make_sym

def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons

plot_dir = '/scratch/rg419/plots/itcz_plots/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

def isca_monsoon_hm(run):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')

    mn_dic = month_dic(1)
    tickspace = [1, 19, 37, 55]
    labels = ['1st Jan', '1st Apr', '1st Jul', '1st Oct']     
    
    # Set figure parameters
    rcParams['figure.figsize'] = 15, 6
    rcParams['font.size'] = 14
    
    # Start figure with 8 subplots
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    def precip_plot(data, ax, lons, title, levels=np.arange(0.,21.,2.)):
        lons = pick_lons(data,lons)
        f1 = (data.precipitation*86400.).sel(lon=lons).mean('lon').plot.contourf(ax=ax, x='xofyear', y='lat', levels = levels, add_labels=False, extend='max', cmap='Blues', add_colorbar=False)
        ax.set_ylim(-45,45)
        ax.grid(True,linestyle=':')
        ax.set_yticks(np.arange(-30.,31.,30.))
        ax.set_title(title)
        return f1
    
    f1 = precip_plot(data, ax1, [0.,10.], '0-10 E')
    precip_plot(data, ax2, [80.,100.], '80-100 E')
    #precip_plot(data, ax3, [125.,145.], '125-145 E')
    precip_plot(data, ax3, [170.,180.], '170-180 E')
    precip_plot(data, ax4, [180.,190.], '180-190 E')
    #precip_plot(data, ax5, [215.,235.], '215-235 E')
    precip_plot(data, ax5, [260.,280.], '260-280 E')
    #precip_plot(data, ax7, [305.,325.], '305-325 E')
    precip_plot(data, ax6, [350.,0.], '350-0 E')
    
    ax1.set_ylabel('Latitude')
    ax5.set_ylabel('Latitude')
    for ax in axes:
        ax.set_xticks(tickspace)
        ax.set_xticklabels(labels,rotation=25)
    
    plt.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.2)
    cb1=fig.colorbar(f1, ax=axes, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.15, aspect=60, shrink=0.5)
    cb1.set_label('Precipitation, mm/day')
    
    # Save as a pdf
    plt.savefig(plot_dir + 'precip_monsoons_' + run + '_diffsplit.pdf', format='pdf')
    plt.close()
    
    data.close()


def ap_monsoon_hm(run, title):
    #Produce identical plot to above, but averaged over all latitudes for a uniform mld run
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')

    mn_dic = month_dic(1)
    tickspace = [1, 19, 37, 55]
    labels = ['1st Jan', '1st Apr', '1st Jul', '1st Oct']     
    
    # Set figure parameters
    rcParams['figure.figsize'] = 5, 3.5
    rcParams['font.size'] = 14
    
    fig, ax1 = plt.subplots()
    try:
        f1 = (data.precipitation*86400.).mean('lon').plot.contourf(ax=ax1, x='xofyear', y='lat', levels = np.arange(0.,21.,2.), add_labels=False, extend='max', cmap='Blues', add_colorbar=False)
    except:
        f1 = ((data.condensation_rain + data.convection_rain)*86400.).mean('lon').plot.contourf(ax=ax1, x='xofyear', y='lat', levels = np.arange(0.,21.,2.), add_labels=False, extend='max', cmap='Blues', add_colorbar=False)
    ax1.set_ylim(-45,45)
    ax1.grid(True,linestyle=':')
    ax1.set_yticks(np.arange(-30.,31.,30.))
    ax1.set_title(title)
    ax1.set_ylabel('Latitude')
    ax1.set_xticks(tickspace)
    ax1.set_xticklabels(labels,rotation=25)
    
    plt.subplots_adjust(left=0.18, right=0.97, top=0.9, bottom=0.2, hspace=0.3, wspace=0.2)
    #cb1=fig.colorbar(f1, ax=ax1, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.15, aspect=60, shrink=0.5)
    #cb1.set_label('Precipitation, mm/day')
    
    # Save as a pdf
    plt.savefig(plot_dir + 'precip_monsoons_' + run + '.pdf', format='pdf')
    plt.close()
    
    data.close()
    

#ap_monsoon_hm('ap_2', '2m Slab Ocean')
#ap_monsoon_hm('mld_20', '20m Slab Ocean')
isca_monsoon_hm('half_shallow')
#isca_monsoon_hm('half_shallow_5')
#isca_monsoon_hm('half_shallow_10')
#isca_monsoon_hm('half_nh_shallow')
#isca_monsoon_hm('q_shallow')
#isca_monsoon_hm('3q_shallow')



