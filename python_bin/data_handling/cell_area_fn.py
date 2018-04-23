import numpy as np
from netCDF4 import Dataset

def cell_area(t_res,base_dir):
	"""read in grid from approriate file, and return 2D array of grid cell areas"""
	resolution_file = Dataset(base_dir + 'src/extra/python/scripts/gfdl_grid_files/t'+str(t_res)+'.nc', 'r', format='NETCDF3_CLASSIC')

	lons = resolution_file.variables['lon'][:]
	lats = resolution_file.variables['lat'][:]

	lonb = resolution_file.variables['lonb'][:]
	latb = resolution_file.variables['latb'][:]

	nlon=lons.shape[0]
	nlat=lats.shape[0]

	area_array = np.zeros((nlat,nlon))

	for i in np.arange(len(lons)):
	    for j in np.arange(len(lats)):
	    	area_array[j,i] = np.absolute(np.radians(lonb[i+1]-lonb[i])*np.cos(np.radians(lats[j]))*np.radians(latb[j+1]-latb[j]))

	return area_array


if __name__ == "__main__":

	# specify resolution
	t_res = 42
	# specify base dir
	base_dir= '/scratch/rg419/GFDL_model/GFDLmoistModel/'
	#return area_array
	area_array=cell_area(t_res,base_dir)
