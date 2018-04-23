""" Investigate Emanuel 1995 critical subcloud entropy curvatures """

from data_handling import month_dic
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
from physics import gradients as gr
from climatology import precip_centroid

plot_dir = '/scratch/rg419/plots/subcloud_entropy/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)
    
def pick_lons(data, lonin):
    #Find index range covering specified longitudes
    if lonin[1]>lonin[0]:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] and data.lon[i] < lonin[1]]
    else:
        lons = [data.lon[i] for i in range(len(data.lon)) if data.lon[i] >= lonin[0] or data.lon[i] < lonin[1]]
    return lons
    
def subcloud_entropy_gradient(data, lonin=[-1.,361.]):
    
    # declare constants
    cp = 287.04/2*7
    L = 2.500e6
    
    # pick longitudes to look at
    lons = pick_lons(data, lonin)
    
    # Evaluate equivalent potential temperature
    ept = ((data.temp + L/cp*data.sphum)*(1000./data.pfull)**(2./7.)).sel(pfull=850., lon=lons).mean('lon').drop('pfull')    
    data['ept'] = ept
    
    # moist entropy, s = cp ln(ept)
    s = cp*np.log(ept)
    
    # entropy gradient
    dsdy = gr.ddy(s, vector=False)
    
    return dsdy


def grad_plot(run, lonin=[-1.,361.]):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    dsdy = subcloud_entropy_gradient(data, lonin=lonin)

    dsdy.plot.contourf(x='xofyear', y='lat', levels=np.arange(-6.e-5, 6.1e-5, 5e-6))
    try:
        (data.precipitation*86400.).mean('lon').plot.contour(x='xofyear', y='lat', colors='k', levels = np.arange(2.,15.,2.))
    except:
        data['precipitation'] = (data.convection_rain + data.condensation_rain)*86400.
        data.precipitation.mean('lon').plot.contour(x='xofyear', y='lat', colors='k', levels = np.arange(2.,15.,2.))
    
    plt.savefig(plot_dir + 'entropy_grad_' + run + '.pdf', format='pdf')
    plt.close()
    

def grad_line_plot(run, ax_in, lonin=[-1.,361.]):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    precip_centroid(data)
    
    dsdy = 100000.* subcloud_entropy_gradient(data, lonin=lonin)

    dsdy[:,31:33].mean('lat').plot(ax=ax_in, color='k')
    
    #ax = plt.gca()    
    ax_twin = ax_in.twinx()
    
    data.p_cent.plot(ax=ax_twin, color='b')
    
    ax_in.set_xlabel('')
    ax_in.set_ylabel('')
    ax_twin.set_ylabel('')
    ax_twin.set_ylim([-30,30])
    ax_in.set_ylim([-8., 8.])
    ax_in.set_yticks([-8.,-4.,0.,4.,8.])
    ax_twin.spines['right'].set_color('blue')
    plt.tight_layout()
    
    return ax_twin



def pcent_egrad_scatter(run, ax_in, lonin=[-1.,361.]):
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/' + run + '.nc')
    
    precip_centroid(data)
    
    dsdy = 100000.* subcloud_entropy_gradient(data, lonin=lonin)
    
    ax_in.plot(dsdy[:,31:33].mean('lat'), data.p_cent, 'xk')
        
    ax_in.set_xlabel('')
    ax_in.set_ylabel('')
    ax_in.set_ylim([-30,30])
    ax_in.set_xlim([-8,8])
    
    
    

if __name__ == "__main__":
    
    f, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharey='row')
    axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
    runs = ['sn_1.000', 'zs_sst',
            'ap_2', 'ap_20',
            'rt_0.500', 'rt_2.000',
            'sn_0.500', 'sn_2.000']
    
    axes_twin=[]
    for i in range(8):
        ax_twin = grad_line_plot(runs[i], axes[i])
        axes_twin.append(ax_twin)
        
    for i in range(6):
        axes[i].set_xlim([0,72])
    
    for i in range(0,7,2):
        axes_twin[i].set_yticklabels('')
    
    axes[6].set_xlim([0,36])
    axes[7].set_xlim([0,144])
    
    plt.savefig(plot_dir + 'entropy_grad_line.pdf', format='pdf')
    plt.close()
    
    #ax2.set_ylabel('Precipitation centroid')
    #ax.set_xlabel('Pentad')
    #ax.set_ylabel('Subcloud entropy gradient')
    
    
    
    f2, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, sharey='row', sharex='col')
    axes = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
    runs = ['sn_1.000', 'zs_sst',
            'ap_2', 'ap_20',
            'rt_0.500', 'rt_2.000',
            'sn_0.500', 'sn_2.000']
    
    axes_twin=[]
    for i in range(8):
        ax_twin = pcent_egrad_scatter(runs[i], axes[i])
        
        
    plt.savefig(plot_dir + 'entropy_grad_scatter.pdf', format='pdf')
    plt.close()
    
    #ax2.set_ylabel('Precipitation centroid')
    #ax.set_xlabel('Pentad')
    #ax.set_ylabel('Subcloud entropy gradient')
    