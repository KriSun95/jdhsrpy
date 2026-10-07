"""
===================
Using Event Filters
===================

This example shows tusing the event list filters found in ``jdhsrpy.list_filters``.
"""

import ntpath
import os
import urllib.request

import matplotlib.pyplot as plt
from astropy.visualization import time_support

from jdhsrpy import TEST_DATA_LOCATION
from jdhsrpy.list_filters import by_detector, by_energy
from jdhsrpy.nustar_evt import NustarEvt, draw_grid
from jdhsrpy.visualize import time_profile_plot

# %%
# Download an example "sunpos" NuSTAR EVT file.

url = f"{TEST_DATA_LOCATION}/test_nustar_data/nu20619003001A06_cl_sunpos.evt"
filename = os.path.join(os.getcwd(), ntpath.split(url)[1])
if not os.path.isfile(filename):
    urllib.request.urlretrieve(url, filename)

# %%
# Use the class to load in and work with the file.

nustar_object = NustarEvt(evt_filename=filename)

# %%
# All of the image and time profile methods and functions that depend on
# the event list accept an event list as an input. If one isn't given
# then the default cleaned one is used.
#
# If a user wants to filter the data then all they must do is filter the
# event list and pass it.
#
# There are two event lists a user might be interested in. For the vast
# majority of cases ``nustar_object.cleaned_evt_data`` is going to work.
# This even list has already bee filtered to remove bad pixel data,
# only includes grade 0 counts, and to only include counts with energies
# between 2.5 keV and 80 keV.
#
# The original event list directly from the EVT file is stored in
# ``nustar_object.evt_data``.
#
# The ``jdhsrpy.list_filters`` module has a number of event list filters
# and a few will be shown here.

# %%
# Filter by detector
# ------------------
#
# Let's create a time profile from the counts from each NuSTAR detector.

time_support(format="unix_tai")
plt.figure()
for det in range(4):
    timesd, ctd = nustar_object.count_time_profile_array(
        event_data=by_detector(nustar_object.cleaned_evt_data, det)
    )
    axes = time_profile_plot(timesd, ctd, label=f"Det{det}")
plt.title(f"FPM{nustar_object.fpm} time profile - by detector")
plt.xticks(rotation=30, ha="right")
plt.ylabel("Counts")
plt.xlabel("Time")
plt.legend()
plt.show()

# %%
# As previously stated, the image methods can also accept a new event
# list so let's try the ``field_of_view_map`` method and only look at
# detector 0.

# obtain a Sunpy map object of the NuSTAR observation
m = nustar_object.field_of_view_map(
    event_data=by_detector(nustar_object.cleaned_evt_data, 0)
)
# now use it in plotting, this is now just Sunpy API stuff
fig = plt.figure()
ax = fig.add_subplot(projection=m)
m.plot(axes=ax)
draw_grid(m, ax)
plt.show()

# %%
# Filter by energy
# ----------------
#
# Let's create a time profile from the counts over a few energy ranges.

# choose 3 energy ranges in keV
energy_ranges = [[3, 4], [4.5, 6], [6, 10]]

time_support(format="unix_tai")
plt.figure()
for er in energy_ranges:
    timese, cte = nustar_object.count_time_profile_array(
        event_data=by_energy(
            nustar_object.cleaned_evt_data, energy_low=er[0], energy_high=er[1]
        )
    )
    axes = time_profile_plot(timese, cte, label=f"{er[0]}-{er[1]} keV")
plt.title(f"FPM{nustar_object.fpm} energy time profiles")
plt.xticks(rotation=30, ha="right")
plt.ylabel("Counts")
plt.xlabel("Time")
plt.yscale("log")
plt.legend()
plt.show()

# %%
# Let's show off the ``full_disk_map`` map method while filtering for
# counts 3<=E<15 keV.

# obtain a Sunpy map object of the NuSTAR observation
image_range = [3, 15]
m = nustar_object.full_disk_map(
    event_data=by_energy(
        nustar_object.cleaned_evt_data,
        energy_low=image_range[0],
        energy_high=image_range[1],
    )
)
# now use it in plotting, this is now just Sunpy API stuff
fig = plt.figure()
ax = fig.add_subplot(projection=m)
m.plot(axes=ax)
draw_grid(m, ax)
plt.title(f"Energy range: {image_range[0]}-{image_range[1]} keV")
plt.show()
