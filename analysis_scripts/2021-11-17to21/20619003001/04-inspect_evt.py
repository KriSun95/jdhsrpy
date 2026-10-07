import ntpath
import os.path
import sys

from astropy.visualization import time_support

sys.path.insert(
    0, "/Users/kris/Documents/umnPostdoc/projects/analysis/jessie-hsr/jdhsrpy/"
)
import matplotlib.pyplot as plt

from jdhsrpy import list_filters, nustar_evt, screening, utils, visualize

DIRECTORY, FILENAME = ntpath.split(__file__)

CREATE_FILES_PREFLARE = False
CREATE_FILES_FLARE = False

obs_id = "20619003001"


time_1 = "2021-11-20T02:22:10"
time0 = "2021-11-20T02:25:30"
time1 = "2021-11-20T02:28:50"

# plt.figure()
# cf = "/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/hk/nu20619003001_chu123.fits"
# chu_times, chus, labels = nustar_evt.chu_state_array(cf)
# axes = visualize.chu_plot(chu_times, chus, labels)
# plt.title(f'CHU States of NuSTAR on ' + chu_times[0].strftime('%Y/%m/%d')) #get the date in the title
# plt.xlabel('Start Time - ' + chu_times[0].strftime('%H:%M:%S'))
# plt.ylabel('NuSTAR CHUs')
# visualize.vertical_line_of_time(time_1, c="r", axes=axes)
# visualize.vertical_line_of_time(time0, c="r", axes=axes)
# visualize.vertical_line_of_time(time0, c="g", axes=axes, ls=":")
# visualize.vertical_line_of_time(time1, c="g", axes=axes, ls=":")
# fmt = mdates.DateFormatter('%H:%M')
# axes.xaxis.set_major_formatter(fmt)
# plt.xticks(rotation=30)
# plt.tight_layout()
# plt.savefig(os.path.join(DIRECTORY, f"intermediate_work/chu-time-profile.png"), bbox_inches="tight")
# plt.show()

file_dir = "/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/event_cl/"
file_dir_hk = "/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/hk/"

orig_files = [
    os.path.join(file_dir, f"nu{obs_id}A06_cl_grade0.evt"),
    os.path.join(file_dir, f"nu{obs_id}B06_cl_grade0.evt"),
]
orig_files = [
    os.path.join(file_dir, f"nu{obs_id}A06_cl_sunpos.evt"),
    os.path.join(file_dir, f"nu{obs_id}B06_cl_sunpos.evt"),
]
lvt_files = [
    os.path.join(file_dir_hk, f"nu{obs_id}A_fpm.hk"),
    os.path.join(file_dir_hk, f"nu{obs_id}B_fpm.hk"),
]

for f, lvtf in zip(orig_files, lvt_files):
    nu_obj = nustar_evt.NustarEvt(evt_filename=f)
    fig = plt.figure()
    m = nu_obj.field_of_view_map()
    # m = image_filters.gaussian_filter(nu_obj.field_of_view_map())
    # m = image_filters.deconvolve_with_file(nu_obj.field_of_view_map(), "/usr/local/caldb/data/nustar/fpm/bcf/psf/nuA2dpsfen1_20100101v001.fits")
    # m = nu_obj.full_disk_map()
    # ax = plt.subplot(projection=m, frame_on=False)
    ax = fig.add_subplot(projection=m)
    m.plot(axes=ax)
    nustar_evt.draw_grid(m, ax)
    plt.show()
    break
    times, ct = nu_obj.rate_time_profile_array(lvtf)

    time_support(format="unix_tai")
    plt.figure()
    axes = visualize.time_profile_plot(times, ct)
    visualize.vertical_line_of_time(time_1, c="r", axes=axes)
    visualize.vertical_line_of_time(time0, c="r", axes=axes)
    visualize.vertical_line_of_time(time0, c="g", axes=axes, ls=":")
    visualize.vertical_line_of_time(time1, c="g", axes=axes, ls=":")
    plt.title(f"FPM{nu_obj.fpm} time profile")
    plt.xticks(rotation=30, ha="right")
    plt.ylabel(f"{ct.unit:latex}")
    plt.xlabel("Time")
    plt.savefig(
        os.path.join(DIRECTORY, f"intermediate_work/fpm{nu_obj.fpm}-time-profile.png"),
        bbox_inches="tight",
    )
    plt.show()

    time_support(format="unix_tai")
    plt.figure()
    for det in range(4):
        timesd, ctd = nu_obj.count_time_profile_array(
            list_filters.by_detector(nu_obj.cleaned_evt_data, det)
        )
        axes = visualize.time_profile_plot(timesd, ctd, label=f"Det{det}")
    visualize.vertical_line_of_time(time_1, c="r", axes=axes)
    visualize.vertical_line_of_time(time0, c="r", axes=axes)
    visualize.vertical_line_of_time(time0, c="g", axes=axes, ls=":")
    visualize.vertical_line_of_time(time1, c="g", axes=axes, ls=":")
    plt.title(f"FPM{nu_obj.fpm} time profile - by detector")
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("Counts")
    plt.xlabel("Time")
    plt.legend()
    plt.savefig(
        os.path.join(
            DIRECTORY, f"intermediate_work/fpm{nu_obj.fpm}-time-profile-dets.png"
        ),
        bbox_inches="tight",
    )
    plt.show()

    plt.figure()
    axes = visualize.livetime_plot(*nustar_evt.livetime_array(lvtf))
    visualize.vertical_line_of_time(time_1, c="r", axes=axes)
    visualize.vertical_line_of_time(time0, c="r", axes=axes)
    visualize.vertical_line_of_time(time0, c="g", axes=axes, ls=":")
    visualize.vertical_line_of_time(time1, c="g", axes=axes, ls=":")
    plt.title(f"FPM{nu_obj.fpm} livetime profile")
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("Livetime [%]")
    plt.xlabel("Time")
    plt.savefig(
        os.path.join(
            DIRECTORY, f"intermediate_work/fpm{nu_obj.fpm}-livetime-profile.png"
        ),
        bbox_inches="tight",
    )
    plt.show()

    del nu_obj

    if CREATE_FILES_PREFLARE:
        save_dir = f"/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/event_cl/{utils.only_numbers(time_1)}_to_{utils.only_numbers(time0)}"
        os.makedirs(save_dir, exist_ok=True)
        screening.make_gti_file(
            "/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/event_cl/nu20619003001A06_gti.fits",
            save_name=f"{save_dir}/{utils.only_numbers(time_1)[-4:]}_to_{utils.only_numbers(time0)[-4:]}_gti.fits",
            good_time_interval=[time_1, time0],
            overwrite=True,
        )
        screening.time_filtered_evt_file(
            evt_file=f,
            time_range=[time_1, time0],
            save_dir=f"{save_dir}/",
            overwrite=True,
        )
    if CREATE_FILES_FLARE:
        save_dir = f"/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/event_cl/{utils.only_numbers(time0)}_to_{utils.only_numbers(time1)}"
        os.makedirs(save_dir, exist_ok=True)
        screening.make_gti_file(
            "/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/event_cl/nu20619003001A06_gti.fits",
            save_name=f"{save_dir}/{utils.only_numbers(time0)[-4:]}_to_{utils.only_numbers(time1)[-4:]}_gti.fits",
            good_time_interval=[time0, time1],
            overwrite=True,
        )
        screening.time_filtered_evt_file(
            evt_file=f,
            time_range=[time0, time1],
            save_dir=f"{save_dir}/",
            overwrite=True,
        )

print("Finished.")
