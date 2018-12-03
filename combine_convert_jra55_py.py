import os
import subprocess
import numpy as np
import calendar
import pdb

def join_files(files_in, file_name_out):

    subprocess.call('cdo mergetime '+files_in+' '+file_name_out, shell=True)

def add_to_list(file_in, file_out, nc_file_string_in):
   
    if not os.path.isfile(file_out):
        subprocess.call('cdo -f nc copy '+file_in + ' ' + file_out, shell=True)
    
    if not (file_out in nc_file_string_in):
        nc_file_string=nc_file_string_in+' '+file_out
    else:
        nc_file_string = nc_file_string_in
        
    return nc_file_string        

if __name__=="__main__":

    variable = 'sphum'
    file_prefix = 'anl_p25.051_sphf'
#     file_prefix = 'anl_surf125.002_prmsl'
#     file_prefix = 'anl_p25.039_vvel'
    do_remap = False
    do_day_mean = True
    do_month_mean = True 

    avg_or_daily = '6hourly'
    folder = '1958_2016/'+variable+'_'+avg_or_daily

    years = np.arange(1958,2017)
    months = np.arange(12)+1


    nc_file_string = ''
    if avg_or_daily=='monthly':
        for year_val in years:    
            file_in = './'+folder+'/'+file_prefix+'.'+str(year_val)+'010100_'+str(year_val)+'123100'
            file_out = './'+folder+'/'+file_prefix+'_'+str(year_val)+'010100_'+str(year_val)+'123100.nc'
            add_to_list(file_in, file_out, nc_file_string_in)
            
    if avg_or_daily=='6hourly':
        for year_val in years:    
            for month_val in months:   
                end_month_day = calendar.monthrange(year_val, month_val)[1]

                month_string ='%02d' % month_val
                
                file_in = './'+folder+'/'+file_prefix+'.'+str(year_val)+month_string+'0100_'+str(year_val)+month_string+str(end_month_day)+'18'
                file_out = './'+folder+'/'+file_prefix+'_'+str(year_val)+month_string+'0100_'+str(year_val)+month_string+str(end_month_day)+'18.nc'  
                
                if not(os.path.isfile(file_in) or os.path.isfile(file_out)):
                    file_in = './'+folder+'/'+file_prefix+'.'+str(year_val)+'010100_'+str(year_val)+'123118'
                    file_out = './'+folder+'/'+file_prefix+'_'+str(year_val)+'010100_'+str(year_val)+'123118.nc'
                    
                nc_file_string = add_to_list(file_in, file_out, nc_file_string)
       
    nc_file_out='./'+folder+'/atmos_'+avg_or_daily+'_together.nc'
    if not os.path.isfile(nc_file_out):
        print 'joining files'
        join_files(nc_file_string,nc_file_out)

#     pdb.set_trace()

    if do_day_mean:
        if not os.path.isdir('./1958_2016/'+variable+'_daily/'):
            os.mkdir('./1958_2016/'+variable+'_daily/')
        daymean_file_out = './1958_2016/'+variable+'_daily/atmos_daily_together.nc'
        if not os.path.isfile(daymean_file_out):
            print 'doing daily mean'
            subprocess.call('cdo daymean '+nc_file_out+' '+daymean_file_out,shell=True)        
        
    if do_day_mean and do_month_mean:
        if not os.path.isdir('./1958_2016/'+variable+'_monthly/'):
            os.mkdir('./1958_2016/'+variable+'_monthly/')
        monmean_file_out = './1958_2016/'+variable+'_monthly/atmos_monthly_together.nc'
        if not os.path.isfile(monmean_file_out):
            print 'doing monthly mean from daily data'
            subprocess.call('cdo monmean '+daymean_file_out+' '+monmean_file_out,shell=True)      
    elif (not do_day_mean) and do_month_mean:
        if not os.path.isdir('./1958_2016/'+variable+'_monthly/'):
            os.mkdir('./1958_2016/'+variable+'_monthly/')
        monmean_file_out = './1958_2016/'+variable+'_monthly/atmos_monthly_together.nc'
        if not os.path.isfile(monmean_file_out):
            print 'doing monthly mean from 6hourly data'
            subprocess.call('cdo monmean '+nc_file_out+' '+monmean_file_out,shell=True)      
              
        
    if do_remap:
        remap_file_out = nc_file_out[:-3]+'_2_5.nc'
        if not os.path.isfile(remap_file_out):
            print 'doing remapping'
            subprocess.call('cdo remapbil,r144x73 '+nc_file_out+' '+remap_file_out,shell=True)
        if do_day_mean:
            daymean_remap_file_out = daymean_file_out[:-3]+'_2_5.nc'
            if not os.path.isfile(daymean_remap_file_out):
                print 'doing remapping'
                subprocess.call('cdo remapbil,r144x73 '+daymean_file_out+' '+daymean_remap_file_out,shell=True)
        if do_month_mean:
            monmean_remap_file_out = monmean_file_out[:-3]+'_2_5.nc'
            if not os.path.isfile(monmean_remap_file_out):
                print 'doing remapping'
                subprocess.call('cdo remapbil,r144x73 '+monmean_file_out+' '+monmean_remap_file_out,shell=True)                    