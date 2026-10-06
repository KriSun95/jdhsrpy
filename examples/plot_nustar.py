"""
===================
Using ``NuSTAREvt``
===================

This example shows the useulness of the ``~jdhsrpy.nustar_evt.NustarEvt`` class.
"""

import ntpath
import os

import astropy.units as u
from astropy.visualization import time_support
import matplotlib.pyplot as plt
import urllib.request

from jdhsrpy import TEST_DATA_LOCATION
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
# Time profile
# ------------
#
# The first interesting thing might be to plot the counts time profile.
#
# Notice that the values here are unit aware.

times, counts = nustar_object.count_time_profile_array() 

time_support(format='unix_tai')
plt.figure()
time_profile_plot(times, counts)
plt.title(f"FPM{nustar_object.fpm} time profile")
plt.xticks(rotation=30, ha='right')
plt.ylabel(f"{counts.unit:latex}")
plt.xlabel("Time")
plt.show()

# %%
# Image
# -----
#
# The class also makes it easy to produce images from the EVT file.
#
# There are three user methods that can be used, two to save time and 
# the final more general one called ``nustar_map``.
#
# The general one allows a user to specify the bounds for the map.

# obtain a Sunpy map object of the NuSTAR observation
bottom_left = [-1000, -700] << u.arcsec
top_right = [300, 600] << u.arcsec
m = nustar_object.nustar_map(bottom_left, top_right)
# now use it in plotting, this is now just Sunpy API stuff
fig = plt.figure()
ax = fig.add_subplot(projection=m)
m.plot(axes=ax)
draw_grid(m, ax)
plt.show()

# %%
# The other two save time with ``field_of_view_map`` and ``full_disk_map``,
# each returning a map bounded by the field of view and the full solar 
# disk, respectively.