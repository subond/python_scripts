#Produce plots comparable to Bordoni and Schneider 2008, and also plots of ITCZ latitude based on maximum surface temperature
#Includes functions to load in data for this, to locate the 'ITCZ' latitude, and to produce plots 

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from data_handling import load_year_xr, pentad_dic

def load_bs_vars(run_fol,years,lonrange=['all']):
    #set run name
    year = years[0]
    rundata = load_year_xr(run_fol, year, pinterp=False)
    rundata.coords['pentad'] = (rundata.time // 5) - 71
    mngrp = rundata.ucomp.groupby('pentad').mean(('time'))

    #Initialise arrays to load into       
    cnvp =   xr.DataArray(np.zeros((72,64,len(years))), [('pentad', range(1,73) ), ('lat', rundata.lat), ('year', years )])
    cndp =   xr.DataArray(np.zeros((72,64,len(years))), [('pentad', range(1,73)  ), ('lat', rundata.lat), ('year', years )])
    t_surf =   xr.DataArray(np.zeros((72,64,len(years))), [('pentad', range(1,73)  ), ('lat', rundata.lat), ('year', years )])
    for year in years:
        print year
        rundata = load_year_xr(run_fol, year, pinterp=False)
        rundata.coords['pentad'] = (rundata.time // 5) - 71    
        if lonrange[0] == 'all':   
            cnvp[:,:,year-years[0]]   = rundata.convection_rain.groupby('pentad').mean(('time','lon')) 
            cndp[:,:,year-years[0]]   = rundata.condensation_rain.groupby('pentad').mean(('time','lon')) 
            t_surf[:,:,year-years[0]]   = rundata.t_surf.groupby('pentad').mean(('time','lon')) 
        else:
            lonmin =  np.min([k for k, j in enumerate(rundata.lon) if j >= lonrange[0] ]) 
            #print lonmin, rundata.lon[lonmin]
            lonmax =  np.max([k for k, j in enumerate(rundata.lon) if j <= lonrange[1] ])
            #print lonmax, rundata.lon[lonmax]
            cnvp[:,:,year-years[0]]   = rundata.convection_rain[:,:,lonmin:lonmax].groupby('pentad').mean(('time','lon')) 
            cndp[:,:,year-years[0]]   = rundata.condensation_rain[:,:,lonmin:lonmax].groupby('pentad').mean(('time','lon')) 
            t_surf[:,:,year-years[0]]   = rundata.t_surf[:,:,lonmin:lonmax].groupby('pentad').mean(('time','lon')) 
            
    data_out = xr.Dataset({'totp': (cnvp.mean(('year')) + cndp.mean(('year')))*86400., 't_surf': t_surf.mean(('year')) })
                     
    return data_out
    
    

def itcz_lat_fn(inp_fol, years, lonrange=['all']):

    #read in lats
    data = load_bs_vars(inp_fol, years, lonrange)
    lats = data.lat.values

    #locate peak in zonal mean surface temperature
    itcz_lat = xr.DataArray(np.zeros(72), [('pentad', data.pentad )])
    tmax = data.t_surf.max(('lat')) 
    for i in range(0,72):
        itcz_lat[i] = data.lat[[k for k, j in enumerate(data.t_surf[i,:]) if j == tmax[i]]]

    return itcz_lat




def bs_plot(inp_fol, years, lonrange=['all']):

    data = load_bs_vars(inp_fol, years, lonrange)
    pd_dic = pentad_dic(1)
    tickspace = range(10,80,10)
    labels = [pd_dic[k] for k in tickspace]
    data.totp.plot.contourf(x='pentad', y='lat',levels=np.arange(9.,21.,3.))
    plt.ylim((-40,40))
    plt.xlim((1,73))
    plt.xlabel('Pentad')
    plt.xticks(tickspace,labels,rotation=25)
    plt.ylabel('Latitude')
    plt.title(inp_fol)
    if lonrange[0]=='all':
        plt.savefig('/scratch/rg419/plots/'+inp_fol+'/rainfall_plot_bs.png')
    else:
        plt.savefig('/scratch/rg419/plots/'+inp_fol+'/rainfall_plot_bs_'+str(int(lonrange[0])) + '_' + str(int(lonrange[1])) + '.png')
    plt.clf()
    

def itcz_plot(inp_fols,years, lonrange=['all']):
    
    for inp_fol in inp_fols:
        itcz_lat = itcz_lat_fn(inp_fol, years, lonrange)
        itcz_lat.plot.line() #'x',markersize=3)
    
    pd_dic = pentad_dic(1)
    tickspace = range(10,80,10)
    labels = [pd_dic[k] for k in tickspace]  
    plt.ylim((-60,60))
    plt.xlim((1,73))
    plt.xlabel('Pentad')
    plt.ylabel('Peak SST latitude')
    plt.xticks(tickspace,labels,rotation=25)
    plt.legend(inp_fols, loc='lower right')
    if lonrange[0]=='all':
        plt.savefig('/scratch/rg419/plots/itcz_lat.png')
    else:
        plt.savefig('/scratch/rg419/plots/itcz_lat_'+str(int(lonrange[0])) + '_' + str(int(lonrange[1])) + '.png')
    plt.clf()


    
if __name__ == "__main__":

    bs_plot('aquaplanet_1m_diags', range(2,17),lonrange=[60.,180.])
    bs_plot('aquaplanet_10m_diags', range(2,17),lonrange=[60.,180.])
    bs_plot('flat_10m_diags', range(2,17),lonrange=[60.,180.])
    bs_plot('topo_10m_diags', range(2,17),lonrange=[60.,180.])
    
    inp_fols = ['aquaplanet_1m_diags','aquaplanet_10m_diags','flat_10m_diags','topo_10m_diags']
    itcz_plot(inp_fols,range(2,17),lonrange=[60.,180.])

