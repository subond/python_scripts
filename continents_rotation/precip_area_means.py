""" 03/04/2018 Take area averages of Isca precip over:
India: 70-100E, 5-30N
East Asia: 100-140E, 20-40N
Western North Pacific: 100-170E, 5-20N
America: 240-280E, 5-30N
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from climatology import area_mean


if __name__ == "__main__":
    
    rcParams['figure.figsize'] = 5, 3
    rcParams['font.size'] = 14
    
    plot_dir = '/scratch/rg419/plots/climatology/precip_area_means/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    india_ctrl = area_mean('control_qflux', region='India')
    india_050 = area_mean('control_qflux_0.500', region='India')
    india_075 = area_mean('control_qflux_0.750', region='India')
    india_125 = area_mean('control_qflux_1.250', region='India')
    india_150 = area_mean('control_qflux_1.500', region='India')
    india_175 = area_mean('control_qflux_1.750', region='India')
    india_200 = area_mean('control_qflux_2.000', region='India')
    
    easia_ctrl = area_mean('control_qflux', region='EAsia')
    easia_050 = area_mean('control_qflux_0.500', region='EAsia')
    easia_075 = area_mean('control_qflux_0.750', region='EAsia')
    easia_125 = area_mean('control_qflux_1.250', region='EAsia')
    easia_150 = area_mean('control_qflux_1.500', region='EAsia')
    easia_175 = area_mean('control_qflux_1.750', region='EAsia')
    easia_200 = area_mean('control_qflux_2.000', region='EAsia')
    
    
    
    fig = plt.figure()
    ax = plt.subplot(111)
    india_ctrl.plot(color='k', linewidth=2)
    india_050.plot(linewidth=2)
    india_075.plot(linewidth=2)
    india_125.plot(linewidth=2)
    india_150.plot(linewidth=2)
    india_175.plot(linewidth=2)
    india_200.plot(linewidth=2)
    plt.xlabel('Pentad')
    plt.ylabel('Precipitation, mm/day')
    plt.xlim([1,72])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    #ax.legend(['Control', 'No Tibet'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
    plt.subplots_adjust(left=0.15, right=0.75, top=0.95, bottom=0.2)
    plt.savefig(plot_dir + 'precip_india_rot.pdf', format='pdf')
    plt.close()
    
    fig = plt.figure()
    ax = plt.subplot(111)
    easia_ctrl.plot(color='k', linewidth=2)
    easia_050.plot(linewidth=2)
    easia_075.plot(linewidth=2)
    easia_125.plot(linewidth=2)
    easia_150.plot(linewidth=2)
    easia_175.plot(linewidth=2)
    easia_200.plot(linewidth=2)
    plt.xlabel('Pentad')
    plt.ylabel('Precipitation, mm/day')
    plt.xlim([1,72])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    #ax.legend(['Control', 'No Americas'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
    plt.subplots_adjust(left=0.15, right=0.67, top=0.95, bottom=0.2)
    plt.savefig(plot_dir + 'precip_easia_rot.pdf', format='pdf')
    plt.close()
    
    