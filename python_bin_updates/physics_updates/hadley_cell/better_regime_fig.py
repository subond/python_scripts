"""
19/01/2017
A series of functions to help in plotting max overturning at 500 hPa for cross-equatorial cell vs the latitude of max near surface mse, or vs precip centroid, or vs lat at which cell drops below some threshold
Generally the same as the previous iteration of better_regime_fig (now better_regime_fig_old2) but allowing steady state data to be plotted too.
NB the better_regime_fig_ss.py includes neater labelling, so you might want to use this sometimes. Also update where data is stored now this has moved to disca

"""

import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import numpy as np
import xarray as xr
from climatology import peak_mse, precip_centroid
from pylab import rcParams
from hadley_cell import get_edge_psi
import statsmodels.api as sm


lev = 500. 
    
rcParams['figure.figsize'] = 10, 7
rcParams['font.size'] = 20
rcParams['text.usetex'] = True



def set_vars(data, plottype=None, lonin=[-1.,361.], thresh=0., nh=True):
    # Load up variables to plot, depending on plot type
    edge_loc, psi_max, psi_max_loc = get_edge_psi(data, thresh=thresh, lev=lev, lonin=lonin, nh=nh)
        
    if plottype == 'mse':
        data_mse = peak_mse(data, lonin=lonin)
        vars_out = (data_mse.mse_max_loc, psi_max)
        
    elif plottype == 'pcent':
        data = precip_centroid(data,lonin=lonin)
        vars_out = (data.p_cent, psi_max)

    elif plottype == None:
        vars_out = (edge_loc, psi_max)
    
    return vars_out
        


def load_vars(runs, plottype=None, lonin=None, thresh=0., nh=True):
    # load variables for multiple runs - useful for steady state runs
    vars_out = []
    if lonin==None:
        lonin=[-1.,361.]*len(runs)
    i=0
    for run in runs:
        data = xr.open_dataset('/disca/share/rg419/Data_moist/climatologies/' + run + '.nc')
        run_vars = set_vars(data, plottype=plottype, lonin=lonin[i], thresh=thresh, nh=nh)
        vars_out.append(run_vars)
        i=i+1
    return vars_out



def fit_power_law(edge_loc, psi_max, maxmin, show_coeff=False):
    # edge_loc - latitude of the ITCZ
    # psi_max - max overturning strength
    # maxmin - range of latitudes to use
    # show_coeff - option to print fitting coeffs to screen
    
    # Get latitudes where edge_loc is in the desired range
    lat = (edge_loc <= maxmin[1]) & (edge_loc > maxmin[0])
    edge = np.log( edge_loc[lat] )
    psi_max = np.log( psi_max[lat] )
    
    A = np.array([ edge, np.ones(edge.shape) ])
    
    # Regress psi max as edge_loc + const
    model = sm.OLS(psi_max.values, A.T)
    result=model.fit()
    consts = result.params
    std_err = result.bse
    
    #print result.summary()
    
    if show_coeff:
        # Print the coefficients and standard errors for the fit
        print('=== Coeffs ===')
        print((np.exp(consts[1]), consts[0]))
        print('=== Std Errs ===')
        print((2*std_err[1]*np.exp(consts[1]), 2*std_err[0]))
    
    # Use constants to get the line to plot, return line and constants
    line = np.exp(consts[1]) * np.arange(maxmin[0],maxmin[1]+1)**(consts[0])
    
    return line, consts
    
    
    

def plot_regime(vars_in, varcolor='k', symbol='x', guide=10, do_linefit_l=False, do_linefit_h=False, include_withdrawal=False, show_coeff=False, flipdata=False):
    # varcolor and symbol specify color and symbol for plotting
    # guide helps in locaton of outward and return branch, selected for shem=True
    # do_linefit_l and h, if True fit a power law to the low or high latitude part of the data
    # include_withdrawal, if True plot the whole time series including withdrawal phase
    
    psi_max = vars_in[1]
    edge_loc = vars_in[0]
    
    if flipdata:
        edge_loc = edge_loc*-1.
    
    data_len = edge_loc.shape[0]
    
    # Find the time when the cell edge is furthest south to use as start of time period to plot
    try:
        #i = int(edge_loc[0:50].argmin('xofyear'))
        i = int(edge_loc[0:data_len/2].argmin('xofyear'))
    except:
        i = 0
    # Find when the cell extent is largest and when the cell strength is largest
    # Use the minimum time as a later bound for plotting    
    j1 = int(edge_loc[guide:].argmax('xofyear'))+guide
    j2 = int(psi_max[guide:].argmax('xofyear'))+guide
    if j1<i:  # if edge_loc data is missing j1 may be less than i, so in this case use j2.
        j=j2
    else:
        #j = min([j1,j2])
        j=j2
    print((i,j))
    
    # Plot ITCZ lat on x axis, cell strength on y over this time period
    plt.plot(edge_loc[i:j], psi_max[i:j], symbol, color=varcolor, ms=10, mew=2)
    if include_withdrawal:
        plt.plot(edge_loc[j:], psi_max[j:], symbol, color=varcolor, ms=8, mew=2, alpha=0.5)
        plt.plot(edge_loc[:i], psi_max[:i], symbol, color=varcolor, ms=8, mew=2, alpha=0.5)
    
    # Fit a power law to a given latitude range of the data, and add the line to the plot
    if do_linefit_l:
        line, consts = fit_power_law(edge_loc[i:j], psi_max[i:j], maxmin=[0.5,7.], show_coeff=show_coeff)
        line = np.exp(consts[1]) * np.arange(0.5,7.5,0.5)**(consts[0])
        plt.plot(np.arange(0.5,7.5,0.5), line, varcolor+':', linewidth=2)
    if do_linefit_h:
        line, consts = fit_power_law(edge_loc[i:j], psi_max[i:j], maxmin=[7.,30.], show_coeff=show_coeff)
        plt.plot(np.arange(7.,31.), line, varcolor+':', linewidth=2)




def plot_multiple(var_list, plotname, plot_dir='/scratch/rg419/plots/regime_fig/', logplot=True, g=None, s=None, vc=None, ll=None, lh=None, iw=None, latmax=None, psirange=None, show_coeff=False, flip=None):
    # Overplot regime figs for multiple datasets
    # Flip provides a way to botch things when data is for the wrong half of the year
    
    n = len(var_list)
    
    # If plotting options aren't specified, set to default for all plots
    if s==None:
        s = ['x']*n
    if vc == None:
        vc = ['k']*n
        
    if g==None:
        g = [10]*n
    if ll == None:
        ll = [False]*n
    if lh == None:
        lh = [False]*n
    if iw == None:
        iw = [False]*n
    if flip == None:
        flip = [False]*n
    
    
    for i in range(n):
        try:
            if var_list[i][0].shape==():
                # If only 1 data point exists, this is a steady state run
                if flip[i]:
                    plt.plot(-1.*var_list[i][0], var_list[i][1], color=vc[i], marker = s[i], markersize=10, mew=2.)
                else:
                    plt.plot(var_list[i][0], var_list[i][1], color=vc[i], marker = s[i], markersize=10, mew=2.)
            else:
                # Otherwise need to use plot_regime function to plot it up
                plot_regime(var_list[i], varcolor=vc[i], symbol=s[i], guide=g[i], do_linefit_l=ll[i], do_linefit_h=lh[i], include_withdrawal=iw[i], show_coeff=show_coeff, flipdata=flip[i])
        except:
            print(('Failed for item ' + str(i)))
    
    plt.xlabel('Cell edge')
    plt.ylabel('Max 500 hPa Mass Streamfunction')
    plt.grid(True,linestyle=':')
    if logplot:
        plt.xscale('log')
        plt.yscale('log')
        plt.minorticks_off()
        if latmax != None:
            plt.xlim(0.1, int(latmax) + 5)
            xticklist = [i for i in [1,2,5,10,20,30,40,50] if i < latmax]
            plt.xticks(xticklist)
        if psirange != None:
            plt.ylim(int(psirange[0]), int(psirange[1]) + 50)
            plt.yticks(list(range(int(round(psirange[0]/50)*50), int(round(psirange[1]/50)*50), 100)))
        ax=plt.gca()
        ax.get_xaxis().set_major_formatter(tk.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(tk.ScalarFormatter())
        plt.tight_layout()
    else:
        plt.xlim(0., 25.)
    #    if latmax != None:
    #        plt.ylim([0,int(latmax)+5])
    #    if psirange != None:
    #        plt.xlim([0,int(psirange[1])+50])
    plt.savefig(plot_dir + plotname +'.pdf', format='pdf')
    plt.close()





if __name__ == "__main__":
    
    # Note to future self STOP BEING LAZY, DON'T JUST CREATE PLOTS HERE. You use this script so much it really just makes a mess when you keep deleting things and losing your records!
    
    
    # Oops I did it again, I made plots in the function, will fix this later?
    runs = ['half_shallow', 'half_shallow_5', 'half_shallow_10']
    
    vars_out = load_vars(runs, lonin=[[170.,190.]]*3)
    plot_multiple(vars_out, 'half_shallow_expts_eastcoast', latmax=30., psirange=[50, 550], vc=['C0','C1','C2'], iw=[True]*3)
    vars_out = load_vars(runs, lonin=[[350.,10.]]*3)
    plot_multiple(vars_out, 'half_shallow_expts_westcoast', latmax=30., psirange=[50, 550], vc=['C0','C1','C2'], iw=[True]*3)
    vars_out = load_vars(runs, lonin=[[80.,100.]]*3)
    plot_multiple(vars_out, 'half_shallow_expts_land', latmax=30., psirange=[50, 550], vc=['C0','C1','C2'], iw=[True]*3)
    vars_out = load_vars(runs, lonin=[[260.,280.]]*3)
    plot_multiple(vars_out, 'half_shallow_expts_ocean', latmax=30., psirange=[50, 550], vc=['C0','C1','C2'], iw=[True]*3)
    
    runs = ['half_shallow', 'q_shallow', '3q_shallow']

    #overturning_hm('q_shallow', regions=[[350,10], [35,45], [80,100], [215,235]])
    #overturning_hm('3q_shallow', regions=[[350,10], [125,145], [260,280], [305,325]])
    vars_out = load_vars(runs, lonin=[[170.,190.],[80.,100.],[260.,280.]])
    plot_multiple(vars_out, 'part_shallow_expts_eastcoast', latmax=30., psirange=[50, 550], vc=['C0','C1','C2'], iw=[True]*3)
    vars_out = load_vars(runs, lonin=[[350.,10.]]*3)
    plot_multiple(vars_out, 'part_shallow_expts_westcoast', latmax=30., psirange=[50, 550], vc=['C0','C1','C2'], iw=[True]*3)
    vars_out = load_vars(runs, lonin=[[80.,100.],[35.,55.],[125.,145]])
    plot_multiple(vars_out, 'part_shallow_expts_land', latmax=30., psirange=[50, 550], vc=['C0','C1','C2'], iw=[True]*3)
    vars_out = load_vars(runs, lonin=[[260.,280.],[215.,235.],[305.,325.]])
    plot_multiple(vars_out, 'part_shallow_expts_ocean', latmax=30., psirange=[50, 550], vc=['C0','C1','C2'], iw=[True]*3)
    
    #runs_ss = ['sine_sst_10m_ss_'+str(i) for i in range(180,89,-10)]

    #m = len(runs_ss)
    
    #runs_te = ['sine_sst_10m']
    #runs = runs_te + runs_ss
    #s = ['x'] + ['s']*m
    #flip = [False] + [True]*m
    
    #vars_out = load_vars(runs, nh=False)
    #plot_multiple(vars_out, 'sine_sst_and_ss', s=s, psirange=[50, 550], latmax=30, show_coeff=True, flip=flip)
    
    
    
    #runs_ss = ['sn_1.000_ss_'+str(i) for i in range(180,241,10)]

    #m = len(runs_ss)
    
    #runs_te = ['sn_1.000']
    #runs = runs_te + runs_ss
    #s = ['x'] + ['s']*m
    
    #vars_out = load_vars(runs)
    #plot_multiple(vars_out, 'sn_1.000_and_ss', s=s, psirange=[50, 550], latmax=30, show_coeff=True)




