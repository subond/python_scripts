import subprocess


def make_video(filepattern, output):
    
    #command = 'ffmpeg -framerate 5 -y -i ' + filepattern + ' -c:v libx264 -r 6 -pix_fmt yuv420p -vf scale=3200:-2 ' + output
    command = 'ffmpeg  -framerate 5 -y -start_number 30 -i ' + filepattern + ' -vframes 45 -c:v libx264 -r 6 -pix_fmt yuv420p -vf scale=3200:-2 ' + output
    
    subprocess.call([command], shell=True)
    
            
make_video('/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/none/wind_and_slp_zanom_%02d_none.png', 
                  '/scratch/rg419/plots/zonal_asym_runs/gill_development/half_shallow/video/precip_and_slp_anom.mp4')

#make_video('/scratch/rg419/plots/stretching_video/stretching_pentad_%02d.png', 
#                  '/scratch/rg419/plots/stretching_video/stretching_vid.mp4')

#make_video('/scratch/rg419/plots/horiz_adv_video/horiz_adv_pentad_%02d.png', 
#                  '/scratch/rg419/plots/horiz_adv_video/horiz_adv_vid.mp4')
                  
#make_video('/scratch/rg419/plots/precip_video/rain_and_sst_pentad_%02d.png', 
#                  '/scratch/rg419/plots/precip_video/precip_vid.mp4')
                  
#make_video('/scratch/rg419/plots/wind_video/wind_pentad_%02d.png', 
#                  '/scratch/rg419/plots/wind_video/wind_vid.mp4')
                  