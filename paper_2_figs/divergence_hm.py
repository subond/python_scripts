"""
Functions to look at absolute vorticity evolution in relation to precip centroid latitude (19/04/2018)
"""

import xarray as xr
import sh
import numpy as np
import matplotlib.pyplot as plt
from data_handling_updates import gradients as gr, make_sym, cell_area
from climatology import precip_centroid
from pylab import rcParams
import scipy.interpolate as spint
from hadley_cell import mass_streamfunction

def psi_mean_clim(run):
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    
    area = cell_area(42, '/scratch/rg419/Isca/')
    area_xr = xr.DataArray(area, [('lat', data.lat ), ('lon', data.lon)])
    area_xr = area_xr.mean('lon')
    psi_mean = ((psi*area_xr).sum(('lat'))/area_xr.sum(('lat'))).mean('pfull')
    #psi_mean = psi_mean*-1.f

    return psi_mean
    
def div_plot(run, ax, rot_fac=1., lev=150.):
    '''Plot dvordt or 1/vor * dvordt'''
    
    #Load data
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    #Calculate psi to overplot
    psi = mass_streamfunction(data, a=6376.0e3, dp_in=50.)
    psi /= 1.e9
    psi = make_sym(psi, asym=True)
    
    # Make precip symmetric and find the precip centroid
    data['precipitation'] = make_sym(data.precipitation)
    precip_centroid(data)
    
    # Calculate divergence
    u_dx = gr.ddx(data.ucomp)  # dudx
    v_dy = gr.ddy(data.vcomp)  # dvdy
    div = (u_dx + v_dy).sel(pfull=lev)*86400.
    div = make_sym(div)
    div = div.mean('lon')

    f1 = div.plot.contourf(ax=ax, x='xofyear', y='lat', levels=np.arange(-1.2,1.3,0.2), add_colorbar=False, add_labels=False)
    
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-500.,0.,100.), add_labels=False, colors='0.7', linewidths=2, linestyles='--')
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(0.,510.,100.), add_labels=False, colors='0.7', linewidths=2)
    psi.sel(pfull=500).plot.contour(ax=ax, x='xofyear', y='lat', levels=np.arange(-1000.,1010.,1000.), add_labels=False, colors='0.5', linewidths=2)
    data.p_cent.plot.line(ax=ax, color='k', linewidth=2)
    ax.set_ylim([-60,60])
    ax.set_ylabel('Latitude')
    ax.set_xticks([12,24,36,48,60,72])
    ax.set_yticks([-60,-30,0,30,60])
    ax.set_xlabel('Pentad')
    ax.grid(True,linestyle=':')
    
    box = ax.get_position()
    axColor = plt.axes([box.x0 + box.width*0.92, box.y0+box.y0*0.12, 0.015, box.height])
    cb1=fig.colorbar(f1, cax=axColor, orientation = 'vertical')
    
        

def time_interp(data, times_new):
    fun = spint.interp1d(data.xofyear.values, data.values, axis=0, fill_value='extrapolate', kind='quadratic')
    data_new = fun(times_new)
    data_new = xr.DataArray(data_new, coords=[times_new, data.lat.values, data.lon.values], dims=['xofyear','lat','lon'])
    return data_new
    
    
def lat_interp(data, lats_new):
    fun = spint.interp1d(data.lat.values, data.values, axis=-1, fill_value='extrapolate', kind='quadratic')
    data_new = fun(lats_new)
    data_new = xr.DataArray(data_new, coords=[data.xofyear.values, lats_new], dims=['xofyear','lat'])
    return data_new
    

def div_at_pcent(run, rot_fac=1., lev=150.,lat_bound=45., res=0.01, interp=True):
    '''Calculate the (normalised) vorticity tendency at the precipitation centroid'''
    
    data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    # Interpolate in time if wanted
    if interp:
        times_new = np.arange(1., 72.2, 0.2)
        precip_new = time_interp(data.precipitation, times_new)
        u_new = time_interp(data.ucomp.sel(pfull=lev), times_new)
        v_new = time_interp(data.vcomp.sel(pfull=lev), times_new)
    else:
        precip_new = data.precipitation
        u_new = data.ucomp.sel(pfull=lev)
        v_new = data.vcomp.sel(pfull=lev)
        
    # Find precipitation centroid
    precip_new = make_sym(precip_new)
    data = xr.Dataset({'precipitation': precip_new}, coords=precip_new.coords)
    precip_centroid(data, lat_bound=lat_bound, res=res)
    
    # Calculate divergence
    v_dx = gr.ddx(v_new)  # dvdx
    u_dy = gr.ddy(u_new)  # dudy
    div = (gr.ddx(u_new) + gr.ddy(v_new)).mean('lon') * 86400.
    omega = 7.2921150e-5 * rot_fac
    f = 2 * omega * np.sin(data.lat *np.pi/180)
    vor = (-1.*gr.ddy(u_new)*0. + f).mean('lon') * 86400.
    div = make_sym(div)
    vor = make_sym(vor, asym=True)
    div=div*vor
    div=vor
    # Interpolate divergence in latitude to match precipitation centroid lats
    lats = [data.lat[i] for i in range(len(data.lat)) if data.lat[i] >= -lat_bound and data.lat[i] <= lat_bound]    
    lats_new = np.arange(-lat_bound, lat_bound+res, res)
    div = lat_interp(div.sel(lat=lats), lats_new)
    
    # Get and return dvordt and dvordtvor at the precipitation centroid, as well as the precipitation centroid itself
    div_pcent = [float(div[i,:].sel(lat=data.p_cent.values[i]).values) for i in range(len(data.xofyear))]
    div_pcent = xr.DataArray(np.asarray(div_pcent), coords=[div.xofyear.values], dims=['xofyear'])

    return div_pcent, data.p_cent



if __name__ == "__main__":

    #colors=['b','g','k','r','c','m','y']
    #runs = ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']
    #rot_facs=[0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.]
    #abs_vort_at_pcent_plot(runs, rot_facs, colors)
    #abs_vort_at_pcent_plot(runs, rot_facs, colors, filename='abs_vort_dtvor_pcent_nointerp', interp=False)
    
    #colors=['b','g','k','r','c']
    #runs = ['mld_2.5', 'mld_5', 'sn_1.000', 'mld_15', 'mld_20']
    #rot_facs=[1.]*5
    #abs_vort_at_pcent_plot(runs, rot_facs, colors, filename='abs_vort_dtvor_pcent_mld', interp=False)
    
    #colors=['g','k','r','c','m']
    #runs = ['sn_0.500', 'sn_1.000', 'sn_2.000', 'sn_3.000', 'sn_4.000']
    #rot_facs=[1.]*5
    #period_facs=[0.5,1.,2.,3.,4.]
    #abs_vort_at_pcent_plot(runs, rot_facs, colors, filename='abs_vort_dtvor_pcent_sn', interp=False, period_facs=period_facs)
    
    max_div = []
    rots = np.arange(0.5,2.1,0.25)
    
    for run in ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']:
        i=0
        # Set plotting directory
        plot_dir = '/scratch/rg419/plots/paper_2_figs/divergence/'
        mkdir = sh.mkdir.bake('-p')
        mkdir(plot_dir)

        # Set figure parameters
        rcParams['figure.figsize'] = 5.5, 7.
        rcParams['font.size'] = 14

        # Start figure with 4 subplots
        fig, (ax1, ax2) = plt.subplots(2, 1)
    
        div_plot(run, ax=ax1)

        psi_mean = psi_mean_clim(run)
        div_pcent, p_cent = div_at_pcent(run, interp=False, rot_fac = rots[i])
        dpcentdt = gr.ddt(p_cent)*86400.
        #dpcentdt.plot(ax=ax2, color='k', linewidth=2)
        p_cent.plot(ax=ax2, color='k', linewidth=2)
        print(run, div_pcent.max('xofyear').values)
        #a = div_pcent.where(p_cent==np.min(np.abs(p_cent)), drop=True).values
        #print(run, a)
        max_div.append(div_pcent.max('xofyear').values)
        #max_div.append(a)
        #dpsi_meandt = gr.ddt(psi_mean)*86400./5.
        #psi_mean.plot(ax=ax2, color='g', linewidth=2)
        
        ax2_twin = ax2.twinx()
        #st_pcent.plot(ax=ax2_twin, color='b', linewidth=2)
        div_pcent.plot(ax=ax2_twin, color='b', linewidth=2)
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        ax2.set_xlabel('Pentad')
        ax2.set_ylabel('ITCZ latitude')
        ax2_twin.set_ylabel('Divergence at precip. centroid', color='b')
        ax2.set_xlim([0,72])
        ax2.set_xticks(np.arange(0,73,12))
        ax2.set_yticks(np.arange(-30.,30.1,10.))
        #ax2_twin.set_yticks(np.arange(0.3,0.91,0.1))
        ax2_twin.spines['right'].set_color('blue')
        ax2_twin.yaxis.label.set_color('blue')
        [t.set_color('blue') for t in ax2_twin.yaxis.get_ticklabels()]
        ax2.grid(True,linestyle=':')    
        
        plt.subplots_adjust(left=0.15, right=0.8, top=0.97, bottom=0.1, hspace=0.3)
        
        plt.savefig(plot_dir+'div_' + run + '.pdf', format='pdf')
        plt.close()
        
        i=i+1
    
   # max_div = np.asarray(max_div)
  #  plt.plot(rots, max_div)
   # plt.show()
    
    for run in ['rt_0.500', 'rt_0.750', 'sn_1.000', 'rt_1.250', 'rt_1.500', 'rt_1.750', 'rt_2.000']:
        i=0
        # Set plotting directory
        plot_dir = '/scratch/rg419/plots/paper_2_figs/divergence/'
        mkdir = sh.mkdir.bake('-p')
        mkdir(plot_dir)

        # Set figure parameters
        rcParams['figure.figsize'] = 5.5, 4.
        rcParams['font.size'] = 14


        psi_mean = psi_mean_clim(run)
        div_pcent, p_cent = div_at_pcent(run, interp=False, rot_fac = rots[i])
        dpcentdt = gr.ddt(p_cent)*86400.
        plt.plot(div_pcent, dpcentdt)
    
    plt.show()