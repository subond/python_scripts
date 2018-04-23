import xarray as xarray
import matplotlib.pyplot as plt
import time

start_time = time.time()


# open monthly files and treat as one
wp1 = xarray.open_mfdataset(	['/scratch/sit204/Data_2013/warmpool_cool_1/run%d/atmos_daily.nc' % m for m in range(13, 61)],
	decode_times=False,  # no calendar so tell netcdf lib
	# choose how data will be broken down into manageable
	# chunks.
	chunks={'time': 30,
			'lon': 128//4,
			'lat': 64//2})

wp3 = xarray.open_mfdataset(	['/scratch/sit204/Data_2013/warmpool_cool_2/run%d/atmos_daily.nc' % m for m in range(13, 61)],
	decode_times=False,  # no calendar so tell netcdf lib
	# choose how data will be broken down into manageable
	# chunks.
	chunks={'time': 30,
			'lon': 128//4,
			'lat': 64//2})

read_1_time = time.time()

zmzw1 = wp1.ucomp.mean(('time', 'lon'))

read_2_time = time.time()

zmzw3 = wp3.ucomp.mean(('time', 'lon'))

zmzwdiff=zmzw3-zmzw1

read_3_time = time.time()

data1={'ucomp1':zmzw1, 'ucomp3':zmzw3, 'ucompdiff':zmzwdiff}

read_4_time = time.time()

# contour plot with 25 colour levels
#zmzwdiff.plot.contourf(x='lat', y='pfull', levels=25)
#plt.ylim(wp1.pfull.max(), wp1.pfull.min())
#plt.yscale('log')

#plt.figure()

tsurf1 = wp1.t_surf.mean(('time'))
tsurf3 = wp3.t_surf.mean(('time'))

tsurfdiff=tsurf3-tsurf1

#contour plot with 25 colour levels
tsurfdiff.plot.contourf(x='lon', y='lat', levels=25)

final_time = time.time()

av_list={}

for varz in ['ucomp','vcomp','temp']:

	av_list[varz] = getattr(wp1,varz).mean(('time', 'lon'))



plt.show()
