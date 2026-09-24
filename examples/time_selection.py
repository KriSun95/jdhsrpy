"""
================
Selecting a time
================

This example is in progress.

"""

import numpy as np

#####################################################
#
# Start by creating simulated data and instrument.
# This would all be provided by a given observation.
#
# Can define the photon energies

start, inc = 1.6, 0.04
stop = 80 + inc / 2
ph_energies = np.arange(start, stop, inc)
