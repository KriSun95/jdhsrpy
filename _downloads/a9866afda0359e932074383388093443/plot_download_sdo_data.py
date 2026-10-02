"""
============================================
Downloading SDO Data and Creating Fe18 Proxy
============================================

This example shows a function in the package being used to download SDO
data.

Following this, the example then goes on to create the Fe18 proxy 
described in Del Zanna 2013
"""

import os

import astropy.units as u
import matplotlib.pyplot as plt
from sunpy import config
import sunpy.map
from sunpy.time import parse_time

from jdhsrpy.iron_channel import create_iron18
from jdhsrpy.sdo_download import sdo_download

TIME_FORMAT = config.get("general", "time_format")

# %%
# Download 1 minute's worth of data for the 94 A SDO/AIA channel

file_directory = os.getcwd()
sdo_download(start_time="2018-09-09T10:25:00", 
                end_time="2018-09-09T10:26:00",
                directory=file_directory,
                wave_values=[94, 171, 211],
                get_hmi=False,
                )

# %%
# Can plot one of the files contents using Sunpy

dir_094 = os.path.join(file_directory, "94angstrom")
file_94 = [f for f in os.listdir(dir_094) if f.endswith("fits")]
map_094 = sunpy.map.Map(os.path.join(dir_094, file_94[0]))

plt.figure()
map_094.plot(clip_interval=(1, 99.99)*u.percent)
plt.show()

# %%
# Create the Fe18 files

dir_171 = os.path.join(file_directory, "171angstrom")
dir_211 = os.path.join(file_directory, "211angstrom")
dir_fe18 = os.path.join(file_directory, "fe18")
create_iron18(dir_094, dir_171, dir_211, outdir=dir_fe18, needing_prepped=True)

# %%
# Can plot the result
file_fe18 = [f for f in os.listdir(dir_fe18) if f.endswith("fits")]
map_fe18 = sunpy.map.Map(os.path.join(dir_fe18, file_fe18[0]))
map_fe18.plot_settings.update({"cmap":"Blues_r"})

plt.figure()
map_fe18.plot(clip_interval=(1, 99.99)*u.percent)
plt.title(f"Fe18 {parse_time(map_fe18.date).strftime(TIME_FORMAT)}")
plt.show()
