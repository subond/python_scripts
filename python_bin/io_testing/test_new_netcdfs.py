"""Check updated netcdfs load ok"""

import numpy as np
from data_handling import time_means
import xarray as xr
import matplotlib.pyplot as plt

test = time_means('ap_1_rd', [1,13], filename='atmos_pentad',timeav='pentad')

test.ucomp[71,37,:,:].plot.contourf()
plt.show()