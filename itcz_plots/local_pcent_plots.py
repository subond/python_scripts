''' 
19/07/2018 Plot hms of the precipitation and precip centroid over different longitude ranges
14/08/2018 Modify to allow run, land mask, and averaging regions to be specified. NB should check use of precip_centroid_ll vs precip_centroid
'''
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from climatology import precip_centroid_ll, precip_mse_plot
from data_handling_updates import cell_area_from_xar, gradients as gr, make_sym

plot_dir = '/scratch/rg419/plots/itcz_plots/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

def local_pcent_plots(run, regions=[[350,10], [80,100], [170,190], [260,280]], do_make_sym=True):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    if do_make_sym:
        data['precipitation'] = make_sym(data.precipitation)
    
    data = precip_centroid_ll(data, lat_bound=30.)
    
    def get_lons(lonin, data):
        if lonin[1]>lonin[0]:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
        else:
            lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
        return lons
    
    lons_1 = get_lons(regions[0], data)
    lons_2 = get_lons(regions[1], data)
    lons_3 = get_lons(regions[2], data)
    lons_4 = get_lons(regions[3], data)
    
    dpcentdt_1 = gr.ddt(data.p_cent.sel(lon=lons_1).mean('lon')) * 86400. 
    dpcentdt_2 = gr.ddt(data.p_cent.sel(lon=lons_2).mean('lon')) * 86400. 
    dpcentdt_3 = gr.ddt(data.p_cent.sel(lon=lons_3).mean('lon')) * 86400. 
    dpcentdt_4 = gr.ddt(data.p_cent.sel(lon=lons_4).mean('lon')) * 86400. 
    
    # Set figure parameters
    rcParams['figure.figsize'] = 10, 7
    rcParams['font.size'] = 14
    # Start figure with 4 subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        
    ax1.plot(data.p_cent.sel(lon=lons_1).mean('lon'), dpcentdt_1, 'kx', mew=2, ms=10, alpha=0.5)
    ax1.plot(data.p_cent.sel(lon=lons_1).mean('lon'), dpcentdt_1, alpha=0.3, color='k', linewidth=2)
    ax1.set_title('West coast')
    
    ax2.plot(data.p_cent.sel(lon=lons_2).mean('lon'), dpcentdt_2, 'kx', mew=2, ms=10, alpha=0.5)
    ax2.plot(data.p_cent.sel(lon=lons_2).mean('lon'), dpcentdt_2, alpha=0.3, color='k', linewidth=2)
    ax2.set_title('Land')
    
    ax3.plot(data.p_cent.sel(lon=lons_3).mean('lon'), dpcentdt_3, 'kx', mew=2, ms=10, alpha=0.5)
    ax3.plot(data.p_cent.sel(lon=lons_3).mean('lon'), dpcentdt_3, alpha=0.3, color='k', linewidth=2)
    ax3.set_title('East coast')
    
    ax4.plot(data.p_cent.sel(lon=lons_4).mean('lon'), dpcentdt_4, 'kx', mew=2, ms=10, alpha=0.5)
    ax4.plot(data.p_cent.sel(lon=lons_4).mean('lon'), dpcentdt_4, alpha=0.3, color='k', linewidth=2)
    ax4.set_title('Ocean')
    
    ax3.set_xlabel('ITCZ latitude')
    ax4.set_xlabel('ITCZ latitude')
    ax1.set_ylabel('ITCZ migration rate')
    ax3.set_ylabel('ITCZ migration rate')
    
    for ax in [ax1,ax2,ax3,ax4]:
        ax.grid(True,linestyle=':')
        ax.set_xlim(-25,25)
        ax.set_ylim(-1.0,1.0)
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    # Save as a pdf
    plt.savefig(plot_dir + 'bowties_' + run + '.pdf', format='pdf')
    plt.close()
    
    
    # Start figure with 4 subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    
    precip_mse_plot(data, ax1, lonin=regions[0], lat_bound=30.)
    precip_mse_plot(data, ax2, lonin=regions[1], lat_bound=30.)
    precip_mse_plot(data, ax3, lonin=regions[2], lat_bound=30.)
    precip_mse_plot(data, ax4, lonin=regions[3], lat_bound=30.)
    
    ax1.set_title('West coast')
    ax2.set_title('Land')
    ax3.set_title('East coast')
    ax4.set_title('Ocean')
    
    ax3.set_xlabel('Pentad')
    ax4.set_xlabel('Pentad')
    ax1.set_ylabel('Latitude')
    ax3.set_ylabel('Latitude')
    
    for ax in [ax1,ax2,ax3,ax4]:
        ax.grid(True,linestyle=':')
    
    plt.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1, hspace=0.3, wspace=0.3)
    
    # Save as a pdf
    plt.savefig(plot_dir + 'precip_' + run + '.pdf', format='pdf')
    plt.close()
    
    data.close()



local_pcent_plots('half_nh_shallow', do_make_sym=False)
#local_pcent_plots('q_shallow', regions=[[350,10], [35,45], [80,100], [215,235]])
#local_pcent_plots('3q_shallow', regions=[[350,10], [125,145], [260,280], [305,325]])
#local_pcent_plots('half_shallow_5')
#local_pcent_plots('half_shallow_10')

