"""
====================
Downloading SDO Data
====================

This example shows a function in the package being used to download SDO
data.
"""

import os

import matplotlib.pyplot as plt
import sunpy.map

from jdhsrpy.sdo_download import sdo_download

# %%
# Download a minute's worth of data for the 94 A SDO/AIA channel:

file_directory = os.getcwd()
sdo_download(start_time="2018-09-09T10:25:00", 
                end_time="2018-09-09T10:26:00",
                directory=file_directory,
                wave_values=[94],
                get_hmi=False,
                )

# %%
# Can plot one of the files contents using Sunpy:

file_94 = [f for f in file_directory if f.endswith("fits")]
map_94 = sunpy.map.Map(os.path.join(file_directory, file_94[0]))

plt.figure()
map_94.plot()
plt.show()
