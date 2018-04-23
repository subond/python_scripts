"""Plot the overturning vs cell edge figure for the seasonal cycle runs 08/02/2018"""

from hadley_cell import load_vars, plot_multiple
import matplotlib.pyplot as plt

    
runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
        'sn_2.000', 'sn_4.000', 'sn_8.000' ]


    
vars_cell_edge = load_vars(runs)
vars_pcent = load_vars(runs, plottype='pcent')

for i in range(6):
    plt.plot(vars_pcent[i][0], vars_cell_edge[i][0])
    plt.show()
