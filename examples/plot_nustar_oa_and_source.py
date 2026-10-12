"""
=======================================================
NuSTAR Optical Axis, Source Location, and Deconvolution
=======================================================

This example shows how to obtain the optical axis information from the
NuSTAR files.

This example will also show how to use the optical axis and source 
location to perform deconvolution correctly on a NuSTAR image.
"""

import ntpath
import os
import urllib.request

import astropy.units as u
import matplotlib.pyplot as plt
from astropy.visualization import time_support

from jdhsrpy import TEST_DATA_LOCATION
from jdhsrpy.nustar_core import NustarEvt, draw_grid

from jdhsrpy.visualize import time_profile_plot

# %%
# Download an example "sunpos" NuSTAR EVT file.

url = f"{TEST_DATA_LOCATION}/test_nustar_data/nu20619003001A06_cl_sunpos.evt"
filename = os.path.join(os.getcwd(), ntpath.split(url)[1])
if not os.path.isfile(filename):
    urllib.request.urlretrieve(url, filename)

url = f"{TEST_DATA_LOCATION}/test_nustar_data/nu20619003001A_oa.fits"
filename_oa = os.path.join(os.getcwd(), ntpath.split(url)[1])
if not os.path.isfile(filename_oa):
    urllib.request.urlretrieve(url, filename_oa)

# %%
# Use the class to load in and work with the file.

nustar_object = NustarEvt(evt_filename=filename)
nustar_object.optical_axis_array(filename_oa)
exit()

# %%
# Time profile
# ------------
#
# The first interesting thing might be to plot the counts time profile.
#
# Notice that the values here are unit aware.

times, counts = nustar_object.count_time_profile_array()

time_support(format="unix_tai")
plt.figure()
time_profile_plot(times, counts)
plt.title(f"FPM{nustar_object.fpm} time profile")
plt.xticks(rotation=30, ha="right")
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
