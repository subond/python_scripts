"""
Load in the vorticity budget terms and produce a lat-time plot at 150 hPa

"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sh
from physics import gradients as gr, model_constants as mc
from pylab import rcParams

def vort_eq_ss(run, lonin=[-1.,361.]):
    
    rcParams['figure.figsize'] = 10, 7.5
    rcParams['font.size'] = 18
    rcParams['text.usetex'] = True
    
    plot_dir = '/scratch/rg419/plots/steady_state_runs/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    #Load in vorticity budget term means
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/vort_eq_'+run+'.nc')    
    print 'vorticity budget data loaded'
    
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]

    # Calculate vertical component of absolute vorticity = f + dv/dx - du/dy
    f = 2 * mc.omega * np.sin(data.lat *np.pi/180)
    v_dx = gr.ddx(data.vcomp)  # dvdx
    u_dy = gr.ddy(data.ucomp)  # dudy
    vor = v_dx - u_dy + f
    
    abs_vort = vor.sel(lon=lons).mean('lon')*86400.
    
    dvordx = gr.ddx(vor)
    dvordy = gr.ddy(vor, vector=False)
    
    horiz_md_mean = -86400.**2. * (data.ucomp * dvordx + data.vcomp * dvordy)
    
    div = gr.ddx(data.ucomp) + gr.ddy(data.vcomp)
    stretching_mean = -86400.**2. * vor * div
    
    transient = (data.horiz_md*86400.**2 + data.stretching*86400.**2 - horiz_md_mean - stretching_mean).values
    
    data['transient'] = (('pfull','lat','lon'), transient )	
    
    horiz_md_hm = horiz_md_mean.sel(lon=lons).mean('lon')
    stretching_hm = stretching_mean.sel(lon=lons).mean('lon')
    transient_hm = data.transient.sel(lon=lons).mean('lon')
    total_hm = data.total.sel(lon=lons).mean('lon')*86400.**2.
    
    
    print 'starting plotting'
    levels = np.arange(-1.5,1.6,0.25)

    # Four subplots
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    plt.set_cmap('RdBu_r')
    #First plot
    f2=horiz_md_hm.plot.contourf(x='lat', y='pfull', levels=levels, ax=ax1, extend = 'both', add_labels=False, add_colorbar=False, yincrease=False)
    ax1.contour(data.lat, data.pfull, abs_vort, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax1.set_xlabel('Latitude')
    ax1.set_xlim(-60,60)
    ax1.set_xticks(np.arange(-60,61,30))
    ax1.grid(True,linestyle=':')

    #Second plot
    stretching_hm.plot.contourf(x='lat', y='pfull', levels=levels, ax=ax2, extend = 'both', add_labels=False, add_colorbar=False, yincrease=False)
    ax2.contour(data.lat, data.pfull, abs_vort, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax2.set_xlim(-60,60)
    ax2.set_xticks(np.arange(-60,61,30))
    ax2.grid(True,linestyle=':')    
    
    #Third plot
    transient_hm.plot.contourf(x='lat', y='pfull', levels=levels, ax=ax3, extend = 'both', add_labels=False, add_colorbar=False, yincrease=False)
    ax3.contour(data.lat, data.pfull, abs_vort, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax3.set_xlabel('Latitude')
    ax3.set_xlim(-60,60)
    ax3.set_xticks(np.arange(-60,61,30))
    ax3.grid(True,linestyle=':')
    
    #Fourth plot
    total_hm.plot.contourf(x='lat', y='pfull', levels = levels, ax=ax4, extend = 'both', add_labels=False, add_colorbar=False, yincrease=False)
    ax4.contour(data.lat, data.pfull, abs_vort, levels=np.arange(-12.,13.,2.), colors='k', linewidths=2, alpha=0.25)
    ax4.set_xlim(-60,60)
    ax4.set_xticks(np.arange(-60,61,30))
    ax4.grid(True,linestyle=':')    
    
    plt.subplots_adjust(right=0.97, left=0.1, top=0.95, bottom=0., hspace=0.2, wspace=0.1)
    #Colorbar
    cb1=f.colorbar(f2, ax=[ax1,ax2,ax3,ax4], use_gridspec=True, orientation = 'horizontal',fraction=0.15, pad=0.1, aspect=30, shrink=0.5)
    cb1.set_label('Vorticity tendency, day$^{-2}$')
        
    if lonin == [-1.,361.]:
        figname = 'vort_budg_ss_' + run + '.pdf'
    else:
        figname = 'vort_budg_ss_' + run + '_' + str(int(lonin[0]))+ '_' + str(int(lonin[1])) + '.pdf'


    plt.savefig(plot_dir + figname, format='pdf')
    plt.close()
    
      
#vort_eq_hm('ss_90.000')
#vort_eq_hm('ss_91.000')
#vort_eq_hm('ss_92.000')
#vort_eq_hm('ss_93.000')
#vort_eq_hm('ss_94.000')
#vort_eq_hm('ss_95.000')
#vort_eq_hm('ss_100.000')
#vort_eq_hm('ss_105.000')




vort_eq_ss('qflux_0_100', [140.,160.])
vort_eq_ss('qflux_10_100', [140.,160.])
vort_eq_ss('qflux_20_100', [140.,160.])
vort_eq_ss('qflux_30_100', [140.,160.])
vort_eq_ss('qflux_0_200', [140.,160.])
vort_eq_ss('qflux_10_200', [140.,160.])
vort_eq_ss('qflux_20_200', [140.,160.])
vort_eq_ss('qflux_30_200', [140.,160.])
vort_eq_ss('qflux_0_300', [140.,160.])
vort_eq_ss('qflux_10_300', [140.,160.])
vort_eq_ss('qflux_20_300', [140.,160.])
vort_eq_ss('qflux_30_300', [140.,160.])