#load in surface temperature in summer and winter and compare. Is the land having any effect?

#import packages
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#set run name
run_fol = 'land_test' #'monsoon_test'


#open files

land_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/exp/' + run_fol + '/input/land.nc'
lfh = Dataset(land_file, 'r', format='NETCDF3_CLASSIC')
land_mask = lfh.variables['land_mask'][:]
lfh.close()

month = ['1','2','3','4','5','6','7','8','9','10','11','12']
t_surf = np.zeros((12,30,64,128))

i=0
for mn in month:
    #nc_file = '/scratch/rg419/Data_moist/' + run_fol + '/np8/run' + mn + '/atmos_daily.nc'
    nc_file = '/scratch/rg419/GFDL_model/GFDLmoistModel/spin_up_restart/' + run_fol + '/np8/run' +mn + '/atmos_daily.nc'
    fh = Dataset(nc_file, mode='r')
    #read in variables
    lons = fh.variables['lon'][:]
    lats = fh.variables['lat'][:]
    pfulls = fh.variables['pfull'][:]
    times = fh.variables['time'][:]
    t_surf[i,:,:,:] = fh.variables['t_surf'][:]
    i=i+1

t_surf_zav = np.repeat(np.reshape(np.mean(t_surf,3),(12,30,64,1)),128,axis=3)
t_surf_zed = t_surf - t_surf_zav

t_surf_zed_mmean = np.squeeze(np.mean(t_surf_zed,1))
t_surf_zed_djf = np.squeeze((t_surf_zed_mmean[0] + t_surf_zed_mmean[1] + t_surf_zed_mmean[11])/3)
t_surf_zed_mam = np.squeeze(np.mean(t_surf_zed_mmean[2:4,:,:],0))
t_surf_zed_jja = np.squeeze(np.mean(t_surf_zed_mmean[5:7,:,:],0))
t_surf_zed_son = np.squeeze(np.mean(t_surf_zed_mmean[8:10,:,:],0))


# Get some parameters for the Stereographic Projection
lon_0 = lons.mean()
lat_0 = lats.mean()

#code from the internet......
#m = Basemap(width=5000000,height=3500000,
#            resolution='l',projection='stere',\
#            lat_ts=40,lat_0=lat_0,lon_0=lon_0)

# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.
lon, lat = np.meshgrid(lons, lats)
m = Basemap(lat_0=lat_0,lon_0=lon_0)
xi, yi = m(lon, lat)



# Plot netcdf Data
def plot_nc_ll(var,tstep,plevel,plotno):
    if tstep < 0:
        varplot = np.squeeze(np.mean(var[:,plevel,:,:],0))
    elif len(var[:].shape) == 4:
        varplot = np.squeeze(var[tstep,plevel,:,:])
    elif len(var[:].shape) == 3:
        varplot = np.squeeze(var[tstep,:,:])
    else:
        print 'what dimension is this...?'

    plt.figure(plotno)
    cs = m.pcolor(xi,yi,varplot)
    plt.contour(xi,yi,land_mask,[0,1],colors='k')
    # Add Grid Lines
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

    # Add Coastlines, States, and Country Boundaries
    #m.drawcoastlines()
    #m.drawstates()
    #m.drawcountries()

    # Add Colorbar
    cbar = m.colorbar(cs, location='bottom', pad="10%")
    cbar.set_label(var.units)

    # Add Title
    plt.title(var.long_name)

    return




# Plot 2d data field
def plot_var_ll(var,tav,plotno):
    if tav < 0:
        varplot = np.squeeze(np.mean(var,0))
    else:
        varplot = var

    plt.figure(plotno)
    cs = m.pcolor(xi,yi,varplot)
    plt.contour(xi,yi,land_mask,[0,1],colors='k')
    # Add Grid Lines
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

    # Add Coastlines, States, and Country Boundaries
    #m.drawcoastlines()
    #m.drawstates()
    #m.drawcountries()

    # Add Colorbar
    cbar = m.colorbar(cs, location='bottom', pad="10%")
    cbar.set_label('Kelvin')

    # Add Title
    plt.title('Surface Temperature')

    return

plot_var_ll(np.squeeze(t_surf_zed[0,:,:,:]),-1,1)
plot_var_ll(np.squeeze(t_surf_zed[5,:,:,:]),-1,2)
#plot_var_ll(temp,1,39,2)
#plot_var_ll(u,29,39,2)
#plot_var_ll(cnv_rain,29,39,3)
#plot_var_ll(slp,29,39,4)

plt.show()


#close files
fh.close()

