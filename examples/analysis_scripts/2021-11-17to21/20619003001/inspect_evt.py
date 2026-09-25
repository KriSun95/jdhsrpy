import os.path
from astropy.visualization import time_support
import sys
from astropy.time import Time
sys.path.insert(0, '/Users/kris/Documents/umnPostdoc/projects/analysis/jessie-hsr/jdhsrpy/')
from jdhsrpy import nustar_evt, visualize, screening, utils
import matplotlib.pyplot as plt

CREATE_FILES = False

obs_id = "20619003001"

file_dir = "/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/event_cl/"

orig_files = [os.path.join(file_dir, f"nu{obs_id}A06_cl_grade0.evt"), 
              os.path.join(file_dir, f"nu{obs_id}B06_cl_grade0.evt")]

for f in orig_files:
    ct, times = nustar_evt.NustarEvt(evt_filename=f).time_profile_array() 
    
    plt.figure()
    plot_times = Time(times, format='unix_tai',scale='utc').datetime
    axes = visualize.time_profile_plot(times, ct)
    time0 = "2021-11-20T02:25:30"
    time1 = "2021-11-20T02:28:50"
    visualize.vertical_line_of_time(time0, c="r", axes=axes)
    visualize.vertical_line_of_time(time1, c="g")
    
    plt.xticks(rotation=30, ha='right')
    plt.show()
    if CREATE_FILES:
        screening.make_gti_file("/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/event_cl/nu20619003001A06_gti.fits", 
                                save_name=f"./{time0}_to_{time1}_gti.fits", 
                                good_time_interval=[time0, time1], 
                                overwrite=True,
                                )
        screening.time_filtered_evt_file(evt_file=f, 
                                        time_range=[time0, time1], 
                                        save_dir="./", 
                                        overwrite=True,
                                        )

print('Finished.')