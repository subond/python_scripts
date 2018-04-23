# Plot the surface temperature and precipitation through the year for the monsoon region

from data_handling_updates import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from climatology import peak_mse, precip_centroid


g = 9.8
cp = 287.04/2*7
L = 2.500e6

def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons
    
def precip_mse_plot(data, ax_in, lonin=[-1.,361.], do_xlabels=False, plot_type=None, precip_contour=8., p_cent=True, mse_max=True, month_labels=True):
    
    lons = pick_lons(data, lonin)
    
    try:
        precip_plot = (data.precipitation*86400.).sel(lon=lons).mean('lon')
        
    except:
        precip_plot = ((data.convection_rain + data.condensation_rain)*86400.).sel(lon=lons).mean('lon')
    
    
    mse_plot = (data.temp*cp + data.sphum*L + data.height*g).mean('lon')/1000.
    
    if plot_type == None:
        # Default case, plot precip with mse overlaid. Choice of whether or not to highlight a specific precip contour
        f1 = precip_plot.plot.contourf(ax=ax_in, x='xofyear', y='lat', levels = np.arange(2.,15.,2.), add_colorbar=False, add_labels=False, extend='max', cmap='Blues')
        if not precip_contour == None:
            precip_plot.plot.contour(ax=ax_in, x='xofyear', y='lat',levels=np.arange(precip_contour-100.,200.,100.), add_labels = False, add_colorbar=False, colors='k', linewidth=2)
        cs = mse_plot.sel(pfull=850.).plot.contour(ax=ax_in, x='xofyear', y='lat', levels=np.arange(200.,401.,10.), add_labels = False, colors='0.7', add_colorbar=False, linewidths=2)
        plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
        if p_cent:
            data = precip_centroid(data,lonin=lonin)
            data.p_cent.plot.line(color='w', ax=ax_in)
            ax_in.set_xlabel('')
        
    elif plot_type == 'precip':
        # No mse, plot precip and precip centroid
        f1 = precip_plot.plot.contourf(ax=ax_in, x='xofyear', y='lat', levels = np.arange(2.,15.,2.), add_colorbar=False, add_labels=False, extend='max', cmap='Blues', linewidth=2)
        if p_cent:
            data = precip_centroid(data,lonin=lonin)
            data.p_cent.plot.line(color='w', ax=ax_in)
            ax_in.set_xlabel('')
            
    elif plot_type == 'mse':
        # Plot mse in colour, overplot max mse
        f1 = mse_plot.sel(pfull=850.).plot.contourf(ax=ax_in, x='xofyear', y='lat', levels=np.arange(200.,401.,10.), add_labels = False, extend='both', add_colorbar=False)
        if mse_max:
            data_mse = peak_mse(data, lonin=lonin)
            data_mse.mse_max_loc.plot.line('k',ax=ax_in)
            ax_in.set_xlabel('')

            
        
    ax_in.set_ylabel('Latitude')
    ax_in.set_ylim(-60,60)
    ax_in.set_yticks(np.arange(-60.,61.,30.))
    ax_in.grid(True,linestyle=':')
    
    if month_labels:
        mn_dic = month_dic(1)
        tickspace = range(13,72,18)
        labels = [mn_dic[(k+5)/6 ] for k in tickspace]
    
        ax_in.set_xlim((1,72))
        ax_in.set_xticks(tickspace)
    
        if do_xlabels:
            ax_in.set_xlabel('')
            ax_in.set_xticklabels(labels,rotation=25)
    
    return f1


if __name__ == "__main__":
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/full_qflux.nc')
    
    plot_dir = '/scratch/rg419/plots/other_monsoons/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #rcParams['figure.figsize'] = 6, 10
    #rcParams['font.size'] = 20
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    f1 = precip_mse_plot(data, ax1, lonin=[345.,45.])
    ax1.set_title('Africa')
    precip_mse_plot(data, ax2, do_xlabels=True, lonin=[240.,270.])
    ax2.set_title('Central America')
    plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0., hspace=0.2)
    #Colorbar
    cb1=fig.colorbar(f1, ax=(ax1, ax2), use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.07, aspect=30)
    cb1.set_label('Precipitation, mm/day')
    
    plt.savefig(plot_dir+'precip_mse_hm.pdf', format='pdf')
    plt.close()        

