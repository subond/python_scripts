"""
Compute t test statistic for the different gradients obtained from the regime fig fits
"""
import numpy as np


def t_test(stats1, stats2):
    
    s = np.sqrt(stats1[1]**2/stats1[2] + stats2[1]**2/stats2[2])
    
    t = (stats2[0] - stats1[0])/s
    
    df = s**4/( (stats1[1]**2/stats1[2])**2/(stats1[2]-1) + (stats2[1]**2/stats2[2])**2/(stats2[2]-1) )
     
    return t, df
    
    

stats1_ap10_500 = [0.09, 0.01, 3.]
stats2_ap10_500 = [0.49, 0.08, 3.]
t, df = t_test(stats1_ap10_500, stats2_ap10_500)
print 'ap10 500', t, df

#stats1_ap10_850 = [0.09, 0.01, 3.]
#stats2_ap10_850 = [0.49, 0.08, 3.]
#t, df = t_test(stats1_ap10_850, stats2_ap10_850)
#print 'ap10 850', t, df

stats1_ap10_400 = [0.05, 0.02, 3.]
stats2_ap10_400 = [0.42, 0.035, 3.]
t, df = t_test(stats1_ap10_400, stats2_ap10_400)
print 'ap10 400', t, df

stats1_ap10_700 = [0.11, 0.02, 3.]
stats2_ap10_700 = [0.24, 0.03, 3.]
t, df = t_test(stats1_ap10_700, stats2_ap10_700)
print 'ap10 700', t, df

stats1_ap10_t90 = [0.14, 0.0005, 3.]
stats2_ap10_t90 = [0.65, 0.105, 3.]
t, df = t_test(stats1_ap10_t90, stats2_ap10_t90)
print 'ap10 t90', t, df

#stats1_ap10_t150 = [0.09, 0.01, 3.]
#stats2_ap10_t150 = [0.49, 0.08, 3.]
#t, df = t_test(stats1_ap10_t150, stats2_ap10_t150)
#print 'ap10 t150', t, df


stats1_full_500 = [0.27, 0.075, 6.]
stats2_full_500 = [0.45, 0.045, 5.]
t, df = t_test(stats1_full_500, stats2_full_500)
print 'full 500', t, df

stats1_full_850 = [0.24, 0.055, 7.]
stats2_full_850 = [0.91, 0.19, 5.]
t, df = t_test(stats1_full_850, stats2_full_850)
print 'full 850', t, df

stats1_full_400 = [0.26, 0.075, 5.]
stats2_full_400 = [0.29, 0.035, 5.]
t, df = t_test(stats1_full_400, stats2_full_400)
print 'full 400', t, df

stats1_full_700 = [0.34, 0.13, 6.]
stats2_full_700 = [0.72, 0.075, 7.]
t, df = t_test(stats1_full_700, stats2_full_700)
print 'full 700', t, df

stats1_full_t90 = [0.24, 0.08, 7.]
stats2_full_t90 = [0.51, 0.055, 6.]
t, df = t_test(stats1_full_t90, stats2_full_t90)
print 'full t90', t, df

stats1_full_t150 = [0.19, 0.265, 3.]
stats2_full_t150 = [0.41, 0.04, 7.]
t, df = t_test(stats1_full_t150, stats2_full_t150)
print 'full t150', t, df


