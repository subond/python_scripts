"""Plot the overturning vs cell edge figure for the seasonal cycle runs 08/02/2018"""

from hadley_cell import load_vars, plot_multiple

    
runs = ['sn_0.250', 'sn_0.500', 'sn_1.000', 
        'sn_2.000', 'sn_4.000', 'sn_8.000' ]

vars_out = load_vars(runs) #, plottype='pcent')
vc = ['m','r','k','y','g','b']
plot_multiple(vars_out, 'seasons',  plot_dir = '/scratch/rg419/plots/paper_2_figs/', vc=vc, psirange=[150, 500], latmax=30, show_coeff=True)


runs = ['rt_0.500', 'rt_0.750', 'sn_1.000', 
        'rt_1.500', 'rt_2.000' ]

vars_out = load_vars(runs) #, plottype='pcent')
vc = ['m','r','k','y','g']
plot_multiple(vars_out, 'rotations',  plot_dir = '/scratch/rg419/plots/paper_2_figs/', vc=vc, psirange=[150, 500], latmax=30, show_coeff=True)


runs = ['ap_2', 'mld_5', 'sn_1.000', 
        'mld_15', 'ap_20' ]

vars_out = load_vars(runs, thresh=0.) #, plottype='pcent')
vc = ['m','r','k','y','g']
plot_multiple(vars_out, 'mlds',  plot_dir = '/scratch/rg419/plots/paper_2_figs/', vc=vc, psirange=[150, 500], latmax=30, show_coeff=True)