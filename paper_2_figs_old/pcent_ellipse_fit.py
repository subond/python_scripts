"""
Find ellipse that fits rate of precip centroid movement. 2/02/2018
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from climatology import precip_centroid
from data_handling_updates import ellipses as el, gradients as gr
from matplotlib.patches import Ellipse

def get_pcent_ellipse(data, check=True):
    # Input: climatology of Isca run precipitation
    # Output: parameters for an ellipse describing the relation between the precipitation centroid and its rate of movement.

    precip_centroid(data)
    dpcentdt = gr.ddt(data.p_cent) * 86400.
    
    lsqe = el.LSqEllipse()
    lsqe.fit([data.p_cent.values, dpcentdt.values])
    center, width, height, phi = lsqe.parameters()
    
    if check:
        plt.close('all')
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(111)
        ax.plot(data.p_cent, dpcentdt, 'ro', label='test data', zorder=1)
        
        ellipse = Ellipse(xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
                   edgecolor='b', fc='None', lw=2, label='Fit', zorder = 2)
        ax.add_patch(ellipse)
        
        plt.legend()
        plt.show()
    
    return center, width, height, phi
    

if __name__ == "__main__":
    
    data = xr.open_dataset('/scratch/rg419/Data_moist/climatologies/sn_0.250.nc')
    
    center, width, height, phi = get_pcent_ellipse(data)
    print center, width, height, phi