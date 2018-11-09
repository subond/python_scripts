'''4/09/2018 Plot onset timing for SCS, Indian, and BoB monsoon onset, for comparison with each other and for validation of code against papers
(JRA, u-850 based metrics)'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sh
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import statsmodels.api as sm

plot_dir = '/scratch/rg419/plots/onset_variability/'
mkdir = sh.mkdir.bake('-p')
mkdir(plot_dir)

# Load in onset dates
onsets_scsm = np.load('/scratch/rg419/python_scripts/onset_variability/jra_onsets_scsm.npy')
onsets_bob  = np.load('/scratch/rg419/python_scripts/onset_variability/jra_onsets_bob_days.npy')
onsets_ism  = np.load('/scratch/rg419/python_scripts/onset_variability/jra_onsets_ism_days.npy')
onsets_scsm = onsets_scsm * 5. - 2.5

# Take mean and standard deviation
scsm_stats = np.mean(np.array(onsets_scsm)), np.std(np.array(onsets_scsm))
bob_stats = np.mean(np.array(onsets_bob)), np.std(np.array(onsets_bob))
ism_stats = np.mean(np.array(onsets_ism)), np.std(np.array(onsets_ism))

years = range(1958,2017)
# Decide whether monsoon is early, normal, or late
def get_timing(onsets, years):
    onset_stats = np.mean(np.array(onsets)), np.std(np.array(onsets))
    print(onset_stats[1])
    enl = []
    enl_no = []
    i=0
    for year in years:
        if onsets[i] > round(onset_stats[0] + 0.75*onset_stats[1]):
            enl.append('Late')
            enl_no.append(2)
        elif onsets[i] < round(onset_stats[0] - 0.75*onset_stats[1]):
            enl.append('Early')
            enl_no.append(0)
        else:
            enl.append('Normal')
            enl_no.append(1)
            
        i=i+1    
    onsets_out = xr.DataArray(onsets, coords={'timing': ('year', enl), 'timing_no': ('year', enl_no), 'year': ('year', years)}, dims=['year'])
    return onsets_out
#print(onsets_scsm.swap_dims({'year': 'timing'}).sel(timing='Normal').year)

onsets_scsm = get_timing(onsets_scsm, years)
onsets_bob = get_timing(onsets_bob, years)
onsets_ism = get_timing(onsets_ism, years)

# Count how often the monsoons are early, normal, or late and make a matrix of these counts
counts = np.zeros((3,3,3))
i=0
for scs in ['Early', 'Normal', 'Late']:
    j=0
    for bob in ['Early', 'Normal', 'Late']:
        k=0
        for ism in ['Early', 'Normal', 'Late']:
            counts[i,j,k] = len(onsets_scsm.year.where( (onsets_scsm['timing']==scs) &
                                          (onsets_ism['timing']==ism) &
                                          (onsets_bob['timing']==bob)
                                           , drop=True))
            k=k+1
        j=j+1
    i=i+1

#print('SCS sum')
#print(sum(counts,0))

#print('BoB sum')
#print(sum(counts,1))

#print('ISM sum')
#print(sum(counts,2))

print(onsets_scsm.year.where( (onsets_ism['timing']=='Early') &
                              #(onsets_ism['timing']=='Normal') #&
                              (onsets_bob['timing']=='Early')
                               , drop=True))

early_bob = onsets_bob.year.where( (onsets_bob['timing']=='Early'), drop=True).values
late_bob = onsets_bob.year.where( (onsets_bob['timing']=='Late'), drop=True).values
early_ism = onsets_ism.year.where( (onsets_ism['timing']=='Early'), drop=True).values
late_ism = onsets_ism.year.where( (onsets_ism['timing']=='Late'), drop=True).values
early_scsm = onsets_scsm.year.where( (onsets_scsm['timing']=='Early'), drop=True).values
late_scsm = onsets_scsm.year.where( (onsets_scsm['timing']=='Late'), drop=True).values

# Make a pop up plot showing how often the different monsoons are early or late
#plt.plot(onsets_scsm.year, onsets_scsm.timing_no,'o')
#plt.plot(onsets_bob.year, onsets_bob.timing_no,'v')
#plt.plot(onsets_ism.year, onsets_ism.timing_no,'x')
#plt.xticks(years)
#plt.yticks([0,1,2])
#plt.grid(True,linestyle=':')
#plt.show()

#Plot all onsets on one axis   
plt.plot(years,onsets_scsm,'o-')
plt.plot(years,onsets_bob,'o-')
plt.plot(years,onsets_ism,'o-')
plt.xlabel('Year')
plt.ylabel('Onset day')
plt.grid(True,linestyle=':')
plt.savefig(plot_dir + 'scs_ism_bob_onsets_jra.pdf', format='pdf')
plt.close()


# Set figure parameters
rcParams['figure.figsize'] = 15, 8
rcParams['font.size'] = 12

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
    
    
# Plot SCSM onset and +-0.75 percentile
ax1.plot(years,onsets_scsm,'o-', color='C0')
ax1.plot([1979,2016],[round(scsm_stats[0])]*2,'-', color='C0')
ax1.fill_between([1979,2016], [round(scsm_stats[0] - 0.75*scsm_stats[1]), round(scsm_stats[0] - 0.75*scsm_stats[1])], 
                              [round(scsm_stats[0] + 0.75*scsm_stats[1]), round(scsm_stats[0] + 0.75*scsm_stats[1])], alpha=0.2, color='C0')
ax1.set_xlabel('Year')
ax1.set_ylabel('Onset day')
ax1.set_title('South China Sea Monsoon Onset')
ax1.text(2005,171,'Mean = %.0f' % scsm_stats[0])
ax1.grid(True,linestyle=':')
#plt.savefig(plot_dir + 'scs_onsets_jra.pdf', format='pdf')
#plt.close()

# Plot BoB onset and +-0.75 percentile
ax2.plot(years,onsets_bob,'o-', color='C1')
ax2.plot([1979,2016],[round(bob_stats[0])]*2,'-', color='C1')
ax2.fill_between([1979,2016], [round(bob_stats[0] - 0.75*bob_stats[1]), round(bob_stats[0] - 0.75*bob_stats[1])], 
                              [round(bob_stats[0] + 0.75*bob_stats[1]), round(bob_stats[0] + 0.75*bob_stats[1])], alpha=0.2, color='C1')
ax2.set_xlabel('Year')
ax2.set_ylabel('Onset day')
ax2.set_title('Bay of Bengal Monsoon Onset')
ax2.text(2005,150,'Mean = %.0f' % bob_stats[0])
ax2.grid(True,linestyle=':')
#plt.savefig(plot_dir + 'bob_onsets_jra.pdf', format='pdf')
#plt.close()

# Plot ISM onset and +-0.75 percentile   
ax3.plot(years,onsets_ism,'o-', color='C2')                          
ax3.plot([1979,2016],[round(ism_stats[0])]*2,'-', color='C2')
ax3.fill_between([1979,2016], [round(ism_stats[0] - 0.75*ism_stats[1]), round(ism_stats[0] - 0.75*ism_stats[1])], 
                              [round(ism_stats[0] + 0.75*ism_stats[1]), round(ism_stats[0] + 0.75*ism_stats[1])], alpha=0.2, color='C2')
ax3.set_xlabel('Year')
ax3.set_ylabel('Onset day')
ax3.set_title('Indian Summer Monsoon Onset')
ax3.text(2005,165.5,'Mean = %.0f' % ism_stats[0])
ax3.grid(True,linestyle=':')
#plt.savefig(plot_dir + 'ism_onsets_jra.pdf', format='pdf')
#plt.close()


# Plot different onsets against one another, along with a best fit line and the 1:1 line

# Bob and SCS
A = np.array([ onsets_bob.values, np.ones(onsets_bob.shape) ])
model = sm.OLS(onsets_scsm.values, A.T)
result=model.fit()
consts = result.params
std_err = result.bse
print('* SCSM onset = (', consts[0], ' +- ', 2*std_err[0], ') * BoB onset + (', consts[1], ' +- ', 2*std_err[1], ') *')
line_bob_scs = consts[1] + np.array([90,155])*(consts[0])

R_bs = np.corrcoef(onsets_bob, onsets_scsm)
R_bi = np.corrcoef(onsets_bob, onsets_ism)
R_si = np.corrcoef(onsets_ism, onsets_scsm)

# SCS and ISM
A = np.array([onsets_scsm.values, np.ones(onsets_scsm.shape) ])
model = sm.OLS(onsets_ism.values, A.T)
result=model.fit()
consts = result.params
std_err = result.bse
print('* ISM onset = (', consts[0], ' +- ', 2*std_err[0], ') * SCSM onset + (', consts[1], ' +- ', 2*std_err[1], ') *')
line_scs_ism = consts[1] + np.array([110,180])*(consts[0])

# Bob and ISM
A = np.array([onsets_bob.values, np.ones(onsets_bob.shape) ])
model = sm.OLS(onsets_ism.values, A.T)
result=model.fit()
consts = result.params
std_err = result.bse
print('* ISM onset = (', consts[0], ' +- ', 2*std_err[0], ') * BoB onset + (', consts[1], ' +- ', 2*std_err[1], ') *')
line_bob_ism = consts[1] + np.array([90,155])*(consts[0])
    
    
ax4.plot(onsets_bob, onsets_scsm, 'xk', ms=10, mew=2)
ax4.plot(onsets_bob.sel(year=early_ism), onsets_scsm.sel(year=early_ism),  'xg', ms=10, mew=2)
ax4.plot(onsets_bob.sel(year=late_ism), onsets_scsm.sel(year=late_ism),  'xr', ms=10, mew=2)
#plt.plot(onsets_bob.sel(year = range(1994,2017)), onsets_scsm.sel(year = range(1994,2017)), 'xr', ms=10, mew=2)
ax4.plot([90,155],line_bob_scs,'-k')
ax4.fill_between([90,155], [round(scsm_stats[0] - 0.75*scsm_stats[1]), round(scsm_stats[0] - 0.75*scsm_stats[1])], 
                              [round(scsm_stats[0] + 0.75*scsm_stats[1]), round(scsm_stats[0] + 0.75*scsm_stats[1])], alpha=0.2, color='C0')
ax4.fill_between([round(bob_stats[0] - 0.75*bob_stats[1]), round(bob_stats[0] + 0.75*bob_stats[1])], [115,115], [175,175], alpha=0.2, color='C1')
ax4.set_ylabel('SCS onset day')
ax4.set_xlabel('BoB onset day')
ax4.set_title('SCSM vs BoB Monsoon Onset')
ax4.plot([0.,180.],[0.,180.],color='0.7', ls='--')
ax4.set_ylim(115,175)
ax4.set_xlim(90,155)
ax4.text(92,171,'R = %.2f' % R_bs[0,1])
#plt.savefig(plot_dir + 'scs_vs_bob_onsets_jra.pdf', format='pdf')
#plt.close()

ax5.plot(onsets_scsm, onsets_ism,  'xk', ms=10, mew=2)
ax5.plot(onsets_scsm.sel(year=early_bob), onsets_ism.sel(year=early_bob),  'xg', ms=10, mew=2)
ax5.plot(onsets_scsm.sel(year=late_bob), onsets_ism.sel(year=late_bob),  'xr', ms=10, mew=2)
#plt.plot(onsets_scsm.sel(year = range(1994,2017)), onsets_ism.sel(year = range(1994,2017)), 'xr', ms=10, mew=2)
ax5.plot([110,180],line_scs_ism,'-k')
ax5.fill_between([115,175], [round(ism_stats[0] - 0.75*ism_stats[1]), round(ism_stats[0] - 0.75*ism_stats[1])], 
                              [round(ism_stats[0] + 0.75*ism_stats[1]), round(ism_stats[0] + 0.75*ism_stats[1])], alpha=0.2, color='C2')
ax5.fill_between([round(scsm_stats[0] - 0.75*scsm_stats[1]), round(scsm_stats[0] + 0.75*scsm_stats[1])], [135,135], [170,170], alpha=0.2, color='C0')
ax5.set_xlabel('SCS onset day')
ax5.set_ylabel('ISM onset day')
ax5.set_title('ISM vs SCSM Monsoon Onset')
ax5.plot([0.,180.],[0.,180.],color='0.7', ls='--')
ax5.set_xlim(115,175)
ax5.set_ylim(135,170)
ax5.text(117,168,'R = %.2f' % R_si[0,1])
#plt.savefig(plot_dir + 'scs_vs_ism_onsets_jra.pdf', format='pdf')
#plt.close()

ax6.plot(onsets_bob, onsets_ism,  'xk', ms=10, mew=2)
ax6.plot(onsets_bob.sel(year=early_scsm), onsets_ism.sel(year=early_scsm),  'xg', ms=10, mew=2)
ax6.plot(onsets_bob.sel(year=late_scsm), onsets_ism.sel(year=late_scsm),  'xr', ms=10, mew=2)
#plt.plot(onsets_bob.sel(year = range(1994,2017)), onsets_ism.sel(year = range(1994,2017)), 'xr', ms=10, mew=2)
ax6.plot([90,155],line_bob_ism,'-k')
ax6.fill_between([90,155], [round(ism_stats[0] - 0.75*ism_stats[1]), round(ism_stats[0] - 0.75*ism_stats[1])], 
                              [round(ism_stats[0] + 0.75*ism_stats[1]), round(ism_stats[0] + 0.75*ism_stats[1])], alpha=0.2, color='C2')
ax6.fill_between([round(bob_stats[0] - 0.75*bob_stats[1]), round(bob_stats[0] + 0.75*bob_stats[1])], [135,135], [170,170], alpha=0.2, color='C1')
ax6.set_xlabel('BoB onset day')
ax6.set_ylabel('ISM onset day')
ax6.set_title('ISM vs BoB Monsoon Onset')
ax6.plot([0.,180.],[0.,180.],color='0.7', ls='--')
ax6.set_xlim(90,155)
ax6.set_ylim(135,170)
ax6.text(92,168,'R = %.2f' % R_bi[0,1])
#plt.savefig(plot_dir + 'ism_vs_bob_onsets_jra.pdf', format='pdf')
#plt.close()

plt.subplots_adjust(left=0.07, right=0.97, top=0.96, bottom=0.07, hspace=0.3, wspace=0.4)

plt.savefig(plot_dir + 'onset_timings_jra.pdf', format='pdf')
plt.close()
