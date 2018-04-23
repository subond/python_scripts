# Plot optical depths as a function of q for: BOG, GEEN 2015, new params

import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams

q = np.arange(0, 0.02, 0.0005)

dtau_bog =  1997.9 * q + 0.8678 
dtau_bog_vic =  1997.9 * q + 0.1627125   

dtau_geen = 351.48*q**0.5 + 0.154925
dtau_geen_win = 1.0814e4*q**2 + 147.11*q + 0.2150

dtau_new = 76.6*q**0.5 + 2.40
dtau_new_win =  2.64e4*q**2 + 215.*q + 0.465


plt.plot(q, dtau_bog, 'b')
plt.plot(q, dtau_bog_vic, 'g')
plt.plot(q, dtau_geen, 'r')
plt.plot(q, dtau_geen_win, 'r--')
plt.plot(q, dtau_new, 'k')
plt.plot(q, dtau_new_win, 'k--')
plt.legend(['BOG', 'BOGV', 'GEEN', 'GEEN_WIN', 'NEW', 'NEW_WIN'], loc=2)
plt.ylabel('dtau/dsigma')
plt.xlabel('Specific humidity, kg/kg')
plt.savefig('/scratch/rg419/plots/radiation_scheme/tau_vs_q')
plt.close