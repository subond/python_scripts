# Plot the streamfunction and jets for the 2m aquaplanet experiment before and after onset, hoping to get a similar result to B+S

#onset pentad:
#topo: 41.25
#flat:41.8
#2m: 37.48
#10m:46.81

from physics import mass_streamfunction
from data_handling import load_year_xr, pentad_dic
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

run_fol = 'aquaplanet_2m'
rundata = load_year_xr('aquaplanet_2m', 21)
years = range(21,41)
pd_dic = pentad_dic(1)

u = xr.DataArray(np.zeros((72,40,64,len(years))), [('pentad', range(1,73) ), ('pfull', rundata.pfull), ('lat', rundata.lat), ('year', years)])
psi = xr.DataArray(np.zeros((64,72,40,len(years))), [('lat', rundata.lat), ('pentad', range(1,73) ), ('pfull', rundata.pfull), ('year', years)])

i=0
for year in years:
    print year
    rundata = load_year_xr(run_fol, year)
    rundata.coords['pentad'] = (rundata.time // 5) -71
    psi_do = mass_streamfunction(rundata, a=6376.0e3)
    u[:,:,:,i] = rundata.ucomp.groupby('pentad').mean(('time','lon'))
    psi[:,:,:,i] = psi_do.groupby('pentad').mean(('time'))
    i=i+1

psi_av = psi.mean('year')/1e9
u_av = u.mean('year')



#data.ucomp.mean(('time','lon')).plot.contourf(x='lat',y='pfull',levels=range(-20,61,5))
cs = psi_av[:,32:36,:].mean('pentad').plot.contour(x='lat',y='pfull', yincrease=False,  colors='k', add_colorbar=False,levels=np.arange(-300.,301.,50.))
u_av[32:36,:,:].mean('pentad').plot.contourf(x='lat',y='pfull',levels=range(-20,61,5), yincrease=False)
plt.xlim(-45,45)
plt.savefig('/scratch/rg419/plots/mom_budg_work/psi_u_before_2m.png')
plt.clf()


cs = psi_av[:,37:41,:].mean('pentad').plot.contour(x='lat',y='pfull', yincrease=False,  colors='k', add_colorbar=False,levels=np.arange(-300.,301.,50.))
u_av[37:41,:,:].mean('pentad').plot.contourf(x='lat',y='pfull',levels=range(-20,61,5), yincrease=False)
plt.xlim(-45,45)
plt.savefig('/scratch/rg419/plots/mom_budg_work/psi_u_after_2m.png')
plt.clf()


g = u_av[26:46,:,:].plot.contourf(x='lat', y='pfull', col='pentad', col_wrap=5,levels=range(-20,61,5), yincrease=False)
plt.xlim(-45,45)
for i, ax in enumerate(g.axes.flat):
    psi_av[:,26+i,:].plot.contour(x='lat', y='pfull', colors='k',ax=ax, levels=np.arange(-300.,301.,50.), add_labels=False, yincrease=False, add_colorbar = False)
    ax.set_title(pd_dic[26+i+1])
 
plt.savefig('/scratch/rg419/plots/mom_budg_work/psi_u_pentads_2m.png')
plt.clf()


