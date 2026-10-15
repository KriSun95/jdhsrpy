"""
===================
Filter By Any Shape
===================

This example shows using more advanced features of the 
``jdhsrpy.list_filters`` region filters.
"""

import ntpath
import os
import urllib.request

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import time_support
from matplotlib import gridspec

from jdhsrpy import TEST_DATA_LOCATION
from jdhsrpy.list_filters import (
    by_general_region,
    by_basic_rectangle_region,
    by_rectangle_region,
)
from jdhsrpy.nustar_core import NustarEvt, draw_grid
from jdhsrpy.utils import arcsec_points_to_shape
from jdhsrpy.visualize import time_profile_plot


# %%
# Download an example "sunpos" NuSTAR EVT file.

url = f"{TEST_DATA_LOCATION}/test_nustar_data/nu20619003001A06_cl_sunpos.evt"
filename = os.path.join(os.getcwd(), ntpath.split(url)[1])
if not os.path.isfile(filename):
    urllib.request.urlretrieve(url, filename)

# %%
# As usual, use the class to load in and work with the file.

nustar_object = NustarEvt(evt_filename=filename)

# %%
# Let's get the default count time profile and image at the start before
# moving on.

bottom_left = [-1000, -500] << u.arcsec
top_right = [100, 425] << u.arcsec

times, ct = nustar_object.count_time_profile_array()
m = nustar_object.nustar_map(bottom_left, top_right)

# %%
# Filter by a rectangle region
# ----------------------------
#
# There are two functions specificly for filtering to a rectangle 
# region. This is probably the most basic and common shape a user might 
# want with which to start.
#
# There is ``by_basic_rectangle_region`` and ``by_rectangle_region`` 
# which should return the exact same result.
# 
# The former is really simple and just compares each counts X and Y 
# position to the maximum and minimum positions of the boxes sides. The 
# latter uses the ``by_general_region`` underneath taking advantage of 
# the fancy `Shapely <https://shapely.readthedocs.io/en/stable/index.html>`__ 
# package.
#
# Here, we will show the use of both region filters.

reg_bottom_left = [-500, 150] << u.arcsec
reg_top_right = [-350, 300] << u.arcsec

basic_rect_event_list = by_basic_rectangle_region(nustar_object.cleaned_evt_data, reg_bottom_left, reg_top_right)
mbr = nustar_object.field_of_view_map(
    event_data=basic_rect_event_list
)
rect_event_list = by_rectangle_region(nustar_object.cleaned_evt_data, reg_bottom_left, reg_top_right)
mr = nustar_object.field_of_view_map(
    event_data=rect_event_list
)

# %% 
# Let's also get the time profile from these region selections.

timesbr, ctbr = nustar_object.count_time_profile_array(
    event_data=basic_rect_event_list
)
timesr, ctr = nustar_object.count_time_profile_array(
    event_data=rect_event_list
)

# %%
#  Also we can use ``~jdhsrpy.list_filters.by_general_region`` filter to 
# produce the same result.

vertices = [reg_bottom_left, 
            [reg_bottom_left[0], reg_top_right[1]]<<u.arcsec, 
            reg_top_right,
            [reg_top_right[0], reg_bottom_left[1]]<<u.arcsec
            ]
polygon = arcsec_points_to_shape(vertices)
general_rect_event_list = by_general_region(nustar_object.cleaned_evt_data, polygon)
mgr = nustar_object.field_of_view_map(
    event_data=general_rect_event_list
)
timesgr, ctgr = nustar_object.count_time_profile_array(
    event_data=general_rect_event_list
)

# %% 
# We can visually see all of these result in the same region and counts 
# are selected in each case.

fig = plt.figure(figsize=(11, 9))
gs = gridspec.GridSpec(2, 2)

axes0 = fig.add_subplot(gs[0, 0], projection=m)
m.plot(axes=axes0, cmap="Greys")
mbr.plot(axes=axes0)
draw_grid(m, axes0)
axes0.set_title(f"``~jdhsrpy.list_filters.by_basic_rectangle_region``")

axes1 = fig.add_subplot(gs[0, 1], projection=m)
m.plot(axes=axes1, cmap="Greys")
mr.plot(axes=axes1)
draw_grid(m, axes1)
axes1.set_title(f"``~jdhsrpy.list_filters.by_rectangle_region``")

axes2 = fig.add_subplot(gs[1, 0], projection=m)
m.plot(axes=axes2, cmap="Greys")
mgr.plot(axes=axes2)
draw_grid(m, axes2)
axes2.set_title(f"``~jdhsrpy.list_filters.by_general_region``")

axes3 = fig.add_subplot(gs[1, 1])
time_profile_plot(times, ct, label="Full observed area", axes=axes3)
time_profile_plot(timesbr, ctbr, label="``~jdhsrpy.list_filters.by_basic_rectangle_region``", ls=":", lw=3, axes=axes3)
time_profile_plot(timesr, ctr, label="``~jdhsrpy.list_filters.by_rectangle_region``", ls="--", lw=2, axes=axes3)
time_profile_plot(timesgr, ctgr, label="``~jdhsrpy.list_filters.by_general_region``", ls="-.", lw=1, axes=axes3)
plt.title(f"FPM{nustar_object.fpm} rectangle time profiles")
plt.xticks(rotation=30, ha="right")
plt.ylabel(f"{ct.unit:latex}")
plt.xlabel("Time")
plt.legend()

plt.tight_layout()
plt.show()

# %%
# It is clear that all the images and time profiles from the region
# selection filters return the same result.
#
# To further make sure all regions are exactly the same we can just 
# subtract the map data and show the sum is ``0``.

diff0 = mr.data-mbr.data
diff1 = mr.data-mgr.data

print(f"``~jdhsrpy.list_filters.by_rectangle_region`` minus ``~jdhsrpy.list_filters.by_basic_rectangle_region`` = {np.sum(diff0)}")
print(f"``~jdhsrpy.list_filters.by_rectangle_region`` minus ``~jdhsrpy.list_filters.by_general_rectangle_region`` = {np.sum(diff1)}")

# %%
# Filter by any shape
# -------------------
#
# The ``~jdhsrpy.list_filters.by_general_rectangle_region`` is very 
# powerful since any region that can be defined by providing vertices of
# a polygon is a region by which we can fitler.
#
# Let's try a few different shapes. 
#
# First, we'll try some sort of a quadrilateral.

vertices0 = [reg_bottom_left, 
             [-550, 320]<<u.arcsec, 
             reg_top_right,
             [-300, 0]<<u.arcsec
             ]
polygon0 = arcsec_points_to_shape(vertices0)
general_rect_event_list0 = by_general_region(nustar_object.cleaned_evt_data, polygon0)
mgr0 = nustar_object.field_of_view_map(
    event_data=general_rect_event_list0
)
timesgr0, ctgr0 = nustar_object.count_time_profile_array(
    event_data=general_rect_event_list0
)

# %%
# Now a triangle.

vertices1 = [[-550, 130]<<u.arcsec, 
             [-400, 300]<<u.arcsec, 
             [-350, 150]<<u.arcsec
             ]
polygon1 = arcsec_points_to_shape(vertices1)
general_rect_event_list1 = by_general_region(nustar_object.cleaned_evt_data, polygon1)
mgr1 = nustar_object.field_of_view_map(
    event_data=general_rect_event_list1
)
timesgr1, ctgr1 = nustar_object.count_time_profile_array(
    event_data=general_rect_event_list1
)

# %%
# Let's be crazy and try a star.

reg_bottom_left = [-500, 150] << u.arcsec
reg_top_right = [-350, 300] << u.arcsec
center = [-425, 225]

vertices2 = [[-475, 125]<<u.arcsec, 
             [-460, 220]<<u.arcsec, 
             [-525, 250]<<u.arcsec,
             [-450, 270]<<u.arcsec,
             [-425, 325]<<u.arcsec, #
             [-400, 270]<<u.arcsec, 
             [-325, 250]<<u.arcsec, 
             [-390, 220]<<u.arcsec,
             [-375, 125]<<u.arcsec,
             [-425, 175]<<u.arcsec,
             ]
polygon2 = arcsec_points_to_shape(vertices2)
general_rect_event_list2 = by_general_region(nustar_object.cleaned_evt_data, polygon2)
mgr2 = nustar_object.field_of_view_map(
    event_data=general_rect_event_list2
)
timesgr2, ctgr2 = nustar_object.count_time_profile_array(
    event_data=general_rect_event_list2
)

# %%
# We can then easily plot the results.

fig = plt.figure(figsize=(11, 9))
gs = gridspec.GridSpec(2, 2)

axes0 = fig.add_subplot(gs[0, 0], projection=m)
m.plot(axes=axes0, cmap="Greys")
mgr0.plot(axes=axes0)
draw_grid(m, axes0)
axes0.set_title("Quadrilateral")

axes1 = fig.add_subplot(gs[0, 1], projection=m)
m.plot(axes=axes1, cmap="Greys")
mgr1.plot(axes=axes1)
draw_grid(m, axes1)
axes1.set_title(f"Triangle")

axes2 = fig.add_subplot(gs[1, 0], projection=m)
m.plot(axes=axes2, cmap="Greys")
mgr2.plot(axes=axes2)
draw_grid(m, axes2)
axes2.set_title(f"Star")

axes3 = fig.add_subplot(gs[1, 1])
time_profile_plot(times, ct, label="Full observed area", axes=axes3)
time_profile_plot(timesgr0, ctgr0, label="Quadrilateral", ls=":", axes=axes3)
time_profile_plot(timesgr1, ctgr1, label="Triangle", ls="--", axes=axes3)
time_profile_plot(timesgr2, ctgr2, label="Star", ls="-.", axes=axes3)
plt.title(f"FPM{nustar_object.fpm} rectangle time profiles")
plt.xticks(rotation=30, ha="right")
plt.ylabel(f"{ct.unit:latex}")
plt.xlabel("Time")
plt.legend()

plt.tight_layout()
plt.show()

# %%
# Therefore, any region that can be defined by its vertices can be used
# to filter the NuSTAR data.
