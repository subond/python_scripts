#Functions for plotting fields up quickly

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

# Do a lat-lon plot
def plot_var_ll(lonsin,lats,var, v=0, plotno=1, pltitle = ' ', cblabel = ' ', land_mask=0, draw_grid=True, draw_coast=False , projin = 'cyl',lon_0=[0.,180.],lat_0=0.,cbarloc='bottom'):

    #define Basemap 
    #lon_0 = lons.mean()
    #lat_0 = lats.mean()
    var,lons = shiftgrid(lon_0[0],var,lonsin)
    lon, lat = np.meshgrid(lons, lats)
    m = Basemap(lat_0=lat_0,lon_0=lon_0[1],projection=projin)
    xi, yi = m(lon, lat)

    #open a specific figure
    plt.figure(plotno)
    #plot contours of variable

    if np.any(v) != 0:
        cs = m.contourf(xi,yi,var,v)
    else:
        cs = m.contourf(xi,yi,var)
    #check if a land-mask was fed in, overplot this if so
    if np.any(land_mask):
        land_mask,lons = shiftgrid(lon_0[0],land_mask,lonsin)
        m.contour(xi,yi,land_mask,[0,1],colors='k')
    # Add Grid Lines
    if draw_grid == True:
        m.drawparallels(np.arange(-60., 61., 30.), labels=[1,0,0,0], fontsize=10)
        m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)

    # Add Coastlines, States, and Country Boundaries
    if draw_coast == True:
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()

    # Add Colorbar
    cbar = m.colorbar(cs, location=cbarloc, pad="10%")
    cbar.set_label(cblabel)

    # Add Title
    plt.title(pltitle)
    #plt.xlabel('Longitude')
    #plt.ylabel('Latitude')
    return



# Do a quiver plot
def plot_quiver_ll(lonsin,lats,var1,var2, plotno=1, pltitle = ' ', cblabel = ' ', land_mask=0, draw_grid=True, draw_coast=False , projin = 'cyl',lon_0=[0.,180.],lat_0=0.):

    #Define Basemap
    #lon_0 = lons.mean()
    #lat_0 = lats.mean()
    var1,lons = shiftgrid(lon_0[0],var1,lonsin)
    var2,lons = shiftgrid(lon_0[0],var2,lonsin)
    lon, lat = np.meshgrid(lons, lats)
    m = Basemap(lat_0=lat_0,lon_0=lon_0[1],projection=projin)
    xi, yi = m(lon, lat)

    #open a specific figure
    plt.figure(plotno)
    # plot quiver of variables
    cs = m.quiver(xi,yi,var1,var2,minshaft=3)
    #check if a land mask was fed in, overplot if so
    if np.any(land_mask):
        land_mask,lons = shiftgrid(lon_0[0],land_mask,lonsin)
        m.contour(xi,yi,land_mask,[0,1],colors='k')
    # Add Grid Lines
    if draw_grid == True:
        m.drawparallels(np.arange(-60., 61., 30.), labels=[1,0,0,0], fontsize=10)
        m.drawmeridians(np.arange(-180., 181., 30.), labels=[0,0,0,1], fontsize=10)

    # Add Coastlines, States, and Country Boundaries
    if draw_coast == True:
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()


    # Add Title
    plt.title(pltitle)
    #plt.xlabel('Longitude')
    #plt.ylabel('Latitude')
    return



# Plot a zonal average
def plot_var_zav(yi,plev,var,v=0, plotno=1, pltitle = ' ', cblabel = ' '):

    #Open up a specific figure
    plt.figure(plotno)
    #If contours are specified, use these, otherwise do automatically
    if np.any(v) != 0:
        cs = plt.contourf(yi,plev,var,v)
    else:
        cs = plt.contourf(yi,plev,var)

    # Add Colorbar
    cbar = plt.colorbar(cs)
    cbar.set_label(cblabel)

    #flip the right way up!
    plt.gca().invert_yaxis()

    # Add Title
    plt.title(pltitle)
    plt.xlabel('Latitude')
    plt.ylabel('Pressure, hPa')

    return



