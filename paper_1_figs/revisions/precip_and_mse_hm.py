# Plot the surface temperature and precipitation through the year for the monsoon region

from data_handling_updates import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh

rcParams['figure.figsize'] = 6, 3.2
rcParams['font.size'] = 17
rcParams['text.usetex'] = True
    
def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons

plot_dir = '/scratch/rg419/plots/paper_1_figs/revisions/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

data = xr.open_dataset('/scratch/rg419/obs_and_reanalysis/gpcp_detrendedpentads.nc')

#lons = pick_lons(data, [60.,150.]) #'Asia'
#lons = pick_lons(data, [15.,60.]) # South Africa
lons = pick_lons(data, [100.,150.]) # East Asia
lons = pick_lons(data, [70.,100.]) # India

data['totp'] = data.precip_clim.sel(lon=lons).mean('lon')

mn_dic = month_dic(1)
tickspace = [13, 31, 49, 68]
labels = ['Mar', 'Jun', 'Sep', 'Dec']    

plevels = np.arange(2.,15.,2.)

# Begin plotting: 3 subplots

#First plot
f1 = data.totp.plot.contourf(x='xofyear', y='lat', levels = plevels, add_labels=False, extend='max', cmap='Blues')
#data.totp.plot.contour(x='xofyear', y='lat',levels=np.arange(-92.,109.,100.), add_labels = False, add_colorbar=False, colors='k')
plt.ylabel('Latitude')
plt.ylim(-60,60)
plt.yticks(np.arange(-60.,61.,30.))
plt.grid(True,linestyle=':')

plt.xlim((1,73))
plt.xlabel('')
plt.xticks(tickspace, labels, rotation=25)
plt.tight_layout() 
plt.savefig(plot_dir+'precip_mse_hm_india.pdf', format='pdf')
plt.close()        

