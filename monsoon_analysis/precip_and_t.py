# Plot the surface temperature and precipitation through the year for the monsoon region

from physics import mombudg_lon_pd_fn
from data_handling import load_year_xr
from pentads import pentad_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


plt.rc('text', usetex=True)
font = {'family' : 'sans-serif','sans-serif':['Helvetica'],
        'weight' : 'bold',
        'size'   : 18}

plt.rc('font', **font)

def plot_weather(inp_fol, years):
    
    year = years[0]
    rundata = load_year_xr(inp_fol, year)
    tot_p = xr.DataArray(np.zeros((72,64,128,len(years))), [('pentad', range(1,73) ), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years)])
    tsurf = xr.DataArray(np.zeros((72,64,128,len(years))), [('pentad', range(1,73) ), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years)])
    #w = xr.DataArray(np.zeros((72,64,128,len(years))), [('pentad', range(1,73) ), ('lat', rundata.lat), ('lon', rundata.lon), ('year', years)])
    
    if not 'aqua' in inp_fol:
        land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/topo_10m/input/land.nc'
        land = xr.open_dataset( land_file)
        land_plot = xr.DataArray(land.land_mask.values, [('lat', rundata.lat), ('lon', rundata.lon)])

    i=0
    for year in years:
        print year
        rundata = load_year_xr(inp_fol, year)
        rundata.coords['pentad'] = (rundata.time // 5) -71
        conv_p = rundata.convection_rain.groupby('pentad').mean(('time'))
        cond_p = rundata.condensation_rain.groupby('pentad').mean(('time'))
        tot_p[:,:,:,i] = (conv_p + cond_p)*86400
        tsurf[:,:,:,i] = rundata.t_surf.groupby('pentad').mean(('time'))
        
        #rundata = load_year_xr(inp_fol, year, pinterp=True)
        #rundata.coords['pentad'] = (rundata.time // 5) -71
        #w[:,:,:,i] = rundata.omega[:,3,:,:].groupby('pentad').mean(('time'))
        
        i=i+1

    tot_p_clim = tot_p.mean(('year'))
    tsurf_clim = tsurf.mean(('year'))
    #w_clim = w.mean(('year'))
    
    pd_dic = pentad_dic(1)
    
    for p in range(0,72):
        print p
        tot_p_clim[p,15:64,21:66].plot.contourf(x='lon', y='lat', add_labels = False, levels=range(0,31,3))
        cs=tsurf_clim[p,15:64,21:66].plot.contour(x='lon', y='lat', colors='w', add_labels = False, add_colorbar=False, levels=range(250,321,5))
        if not 'aqua' in inp_fol:
            land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='w',add_colorbar=False,add_labels=False)
        plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
        plt.ylim(0,90)
        plt.xlim(60,180)
        plt.title(pd_dic[p+1])
        plt.tight_layout()
        plt.savefig('/scratch/rg419/plots/monsoon_analysis/'+inp_fol+'/rain_and_t'+str(p)+'.png')
        plt.close()
        
        #w_clim[p,15:50,21:66].plot.contourf(x='lon', y='lat', add_labels = False, levels=np.arange(-0.3,0.31,0.03))
        #cs=tsurf_clim[p,15:50,21:66].plot.contour(x='lon', y='lat', colors='k', add_labels = False, add_colorbar=False, levels=range(250,321,5))
        #if not 'aqua' in inp_fol:
        #    land_plot.plot.contour(x='lon', y='lat',levels=np.arange(0.,2.,1.), colors='k',add_colorbar=False,add_labels=False)
        #plt.clabel(cs, fontsize=15, inline_spacing=-1, fmt= '%1.0f')
        #plt.ylim(-45,45)
        #plt.xlim(60,180)
        #plt.title(pd_dic[p+1])
        #plt.tight_layout()
        #plt.savefig('/scratch/rg419/plots/monsoon_analysis/'+inp_fol+'/w700_and_t'+str(p)+'.png')
        #plt.close()


#plot_weather('flat_10m',range(11,41))
#plot_weather('topo_10m',range(11,41))
plot_weather('aquamountain_10m',range(11,41))
plot_weather('aquaplanet_2m',range(11,41))
plot_weather('aquaplanet_10m',range(11,41))
