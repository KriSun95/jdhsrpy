"""
=========
REG files
=========

This example shows the useulness of the region file functions in
``~jdhsrpy.screening``.

.. note::

    For reasons, the only regions that can be used with the NuSTAR
    pipeline are circle regions.
"""

import ntpath
import os
import urllib.request

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord

from jdhsrpy import TEST_DATA_LOCATION
from jdhsrpy.io import read_reg_file
from jdhsrpy.nustar_core import NustarEvt, draw_grid, utc_from_nustar_time
from jdhsrpy.shapes import circle_shape
from jdhsrpy.utils import radec_to_helioprojective
from jdhsrpy.visualize import shape_plot

# %%
# Download an example "sunpos" NuSTAR EVT file.

url = f"{TEST_DATA_LOCATION}/test_nustar_data/nu20619003001A06_cl_sunpos.evt"
filename = os.path.join(os.getcwd(), ntpath.split(url)[1])
if not os.path.isfile(filename):
    urllib.request.urlretrieve(url, filename)
url = f"{TEST_DATA_LOCATION}/test_nustar_data/fk5_sexagesimal.reg"
filename_fs = os.path.join(os.getcwd(), ntpath.split(url)[1])
if not os.path.isfile(filename_fs):
    urllib.request.urlretrieve(url, filename_fs)
url = f"{TEST_DATA_LOCATION}/test_nustar_data/icrs_degrees.reg"
filename_id = os.path.join(os.getcwd(), ntpath.split(url)[1])
if not os.path.isfile(filename_id):
    urllib.request.urlretrieve(url, filename_id)

# %%
# Use the class to load in and work with the file.

nustar_object = NustarEvt(evt_filename=filename)

# %%
# Load in a region file
# ---------------------
#
# To load in a region file, say, from SAOImageDS9, dirst a user needs to
# make sure the file was save with either a WCS of FK5 or ICRS.
#
# The loader should be able to handle the values being saved in either
# degrees or sexagesimal regardless of the WCS.

# %%
# FK5 and sexagesimal
# *******************
#
# Let's look at a REG file saved using FK5 and in sexagesimal first.

ra_fs, dec_fs, rad_fs = read_reg_file(filename_fs)

# %%
# Next, let's convert the RA and Dec values to solar X/Y.
#
# For this, we also need a time associate with the conversion. Here,
# I've chosen the mean time of the event list.
#
# Be mindful of the time range of the file used with SAOImageDS9. It
# will likely be the middle of this time range the user wants.

mid_time = np.mean(utc_from_nustar_time(nustar_object.cleaned_evt_data["TIME"]))
sun_x_fs, sun_y_fs = radec_to_helioprojective(ra_fs, dec_fs, mid_time)

# %%
# Now we have everything to define the circle region.

center_fs = (sun_x_fs, sun_y_fs)
radius_fs = rad_fs
circle_region_fs = circle_shape(center_fs, radius_fs)

# %%
# Now we just have to plot the image.
#
# Remember to look at the "Filter By Any Shape" example for different
# uses of the shape objects we use here.

bottom_left = [-1000, -500] << u.arcsec
top_right = [100, 425] << u.arcsec
m = nustar_object.nustar_map(bottom_left, top_right)

# now use it in plotting, this is now just Sunpy API stuff
fig = plt.figure()
ax = fig.add_subplot(projection=m)
m.plot(axes=ax)
shape_plot(circle_region_fs, m.coordinate_frame, color="tab:green")
draw_grid(m, ax)
plt.show()

# %%
# ICRS and degrees
# ****************
#
# We can see that the the WCS ICRS file that had the RA and Dec saved in
# degrees from the same region in SAOImageDS9 give the same circle.

ra_id, dec_id, rad_id = read_reg_file(filename_id)
sun_x_id, sun_y_id = radec_to_helioprojective(ra_id, dec_id, mid_time)
center_id = (sun_x_id, sun_y_id)
radius_id = rad_id
circle_region_id = circle_shape(center_id, radius_id)

fig = plt.figure()
ax = fig.add_subplot(projection=m)
m.plot(axes=ax)
shape_plot(circle_region_fs, m.coordinate_frame, color="tab:green")
shape_plot(circle_region_id, m.coordinate_frame, color="tab:red", ls=":", lw=2)
draw_grid(m, ax)
plt.show()

# %%
# Astropy SkyCoord
# ****************
#
# In order to do all of the coordinates natively in Python. I think the
# following should work.
#
# This has not extensively been tested and does not have heritage like
# the previous approaches that use ``nustar_pysolar`` under the hood.
#
# Use with caution!
#
# The only coordinate system I've been able to get matching the file
# location is the GCRS system.

e_rad = 6371 << u.km  # Earth's radius
dist = 1 << u.au
coord_gcrs = SkyCoord(
    ra=ra_fs,
    dec=dec_fs,
    frame="gcrs",
    distance=dist,
    obstime=mid_time,
    obsgeoloc=(e_rad.value, 0.0, e_rad.value) << e_rad.unit,
    obsgeovel=(0.0, 0.0, 0.0) << (u.m / u.s),
)
gcrs = coord_gcrs.transform_to(m.coordinate_frame)
circle_region_gcrs = circle_shape((gcrs.Tx, gcrs.Ty), radius_fs)

fig = plt.figure()
ax = fig.add_subplot(projection=m)
m.plot(axes=ax)
shape_plot(circle_region_fs, m.coordinate_frame, color="tab:green")
shape_plot(circle_region_gcrs, m.coordinate_frame, color="tab:purple", ls=":", lw=2)
draw_grid(m, ax)
plt.show()

# %%
# Create a region file
# --------------------
