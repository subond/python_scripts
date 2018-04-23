import subprocess


def make_video(filepattern, output):
    
    command = 'ffmpeg -framerate 5 -y -i ' + filepattern + ' -c:v libx264 -r 6 -pix_fmt yuv420p -vf scale=3200:-2 ' + output
    
    subprocess.call([command], shell=True)
    
            
make_video('/scratch/rg419/plots/stretching_video/stretching_pentad_%02d.png', 
                  '/scratch/rg419/plots/stretching_video/stretching_vid.mp4')

make_video('/scratch/rg419/plots/horiz_adv_video/horiz_adv_pentad_%02d.png', 
                  '/scratch/rg419/plots/horiz_adv_video/horiz_adv_vid.mp4')
                  
make_video('/scratch/rg419/plots/precip_video/rain_and_sst_pentad_%02d.png', 
                  '/scratch/rg419/plots/precip_video/precip_vid.mp4')
                  
make_video('/scratch/rg419/plots/wind_video/wind_pentad_%02d.png', 
                  '/scratch/rg419/plots/wind_video/wind_vid.mp4')
                  