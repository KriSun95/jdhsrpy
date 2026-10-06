"""
=======================
NuSTAR Spectral Fitting
=======================

This example shows the spectral fitting of NuSTAR data.
"""

import ntpath
import os
import urllib.request

import matplotlib.pyplot as plt
import numpy as np

from jdhsrpy import TEST_DATA_LOCATION, spectral_fitting
from jdhsrpy.nustar_spectrum import nustar_pha_spectrum
from jdhsrpy.visualize import spectrum_plot

# %%
# Download the spectral files.

url = f"{TEST_DATA_LOCATION}/test_nustar_data/nu20619003001A06_cl_grade0_sr.pha"
filename_pha = os.path.join(os.getcwd(), ntpath.split(url)[1])
if not os.path.isfile(filename_pha):
    urllib.request.urlretrieve(url, filename_pha)

url = f"{TEST_DATA_LOCATION}/test_nustar_data/nu20619003001A06_cl_grade0_sr.arf"
filename_arf = os.path.join(os.getcwd(), ntpath.split(url)[1])
if not os.path.isfile(filename_arf):
    urllib.request.urlretrieve(url, filename_arf)

url = f"{TEST_DATA_LOCATION}/test_nustar_data/nu20619003001A06_cl_grade0_sr.rmf"
filename_rmf = os.path.join(os.getcwd(), ntpath.split(url)[1])
if not os.path.isfile(filename_rmf):
    urllib.request.urlretrieve(url, filename_rmf)

# %%
# This example assumes the user has already obtained the PHA, ARF, and
# RMF files from the NuSTAR data analysis software. See the
# `documentation <https://krisun95.github.io/jdhsrpy/setting_up_nustar_data.html>`__
# for more details.
#
# At any point below, a user can swap to only using `Sunkit-spex <https://sunkit-spex.readthedocs.io/en/latest/>`__
# if the helper functions are too cryptic or annoying.

# %%
# Inspecting the spectrum
# -----------------------
#
# It is possible to inspect the spectrum before continuing with the
# fitting.
#
# This allows a user to get an idea of things like fitting energy ranges
# and which models to use.

energy_bins, flux, flux_err = nustar_pha_spectrum(filename_pha)

plt.figure(figsize=(8, 7))
spectrum_plot(energy_bins, flux, flux_err)
plt.xlim([2, 7])
plt.title("NuSTAR PHA Spectrum")
plt.tight_layout()
plt.show()

# %%
# Setting up the fitting
# ----------------------
#
# We start by loading in the spectral files using ``sunkit-spex`` and
# the helper functions in ``jdhsrpy``.
#
# This example only shows using one NuSTAR observation and only from
# FPMA. To include more, it's just as simple as passing more files in the
# following (e.g. ``spectral_fitting.get_fitter_object(file1, file2, ...)``).

fitter = spectral_fitting.get_fitter_object(filename_pha)

# %%
# The above is the equivalent of the following when using Sunkit-spex
# directly.

from sunkit_spex.legacy.fitting.fitter import Fitter

equivalent_fitter = Fitter(
    pha_file=[
        filename_pha,
    ]
)

# %%
# Choosing the fitting range
# --------------------------
#
# Choose a fitting range in keV.
#
# The fitting range can be changed at any point. E.g., parameters can be
# freed or frozen, the fitting range can be changed, then the fit can be
# run again starting off from where it left off from the previous fit.

fitter.energy_fitting_range = [2.5, 6.5]

# %%
# Choosing the model
# ------------------
#
# There are a lot of models to choose from but let's try a single
# thermal model.
#
# These helper functions are useful as it will check if there are
# multiple spectra being fitted and then will include a constant factor
# automatically in the model if that's the case.

spectral_fitting.set_single_thermal_model(fitter)
print(fitter.params)

# %%
# Setting up the parameters
# -------------------------
#
# For the fitting to start off, we need to provide initial guesses and
# provide bounds for the fitting parameters.

fitter.params["T1_spectrum1"] = {"Value": 3, "Bounds": (1.1, 4)}
fitter.params["EM1_spectrum1"] = {"Value": 9e-1, "Bounds": (5e-1, 5e0)}

# %%
# Fitting the data and plotting
# -----------------------------
#
# The data can now be fitted and plotted.

fitter.fit()

# plot the result
plt.figure(figsize=(8, 7))
axes, res_axes = fitter.plot()
y_max_lim = 1.2 * np.max(fitter.data.loaded_spec_data["spectrum1"]["count_rate"])
axes[0].set_ylim([1e-1, y_max_lim])
axes[0].set_xlim([2, 7])
plt.tight_layout()
plt.show()

# %%
# MCMC analysis
# -------------
#
# TWe can also easily run MCMC analysis.

mcmc_result = fitter.run_mcmc(steps_per_walker=1_000)
fitter.burn_mcmc = 200

# %%
# We can see the log-probability chain of the walkers.

plt.figure()
fitter.plot_log_prob_chain()
plt.ylim([-80, -55])
plt.tight_layout()
plt.show()

# %%
# The corner plot of the walkers.

corner_plot = fitter.corner_mcmc()
plt.tight_layout()
plt.show()

# %%
# And, of course, the final fitted result.

# plot the result
plt.figure(figsize=(8, 7))
axes, res_axes = fitter.plot()
y_max_lim = 1.2 * np.max(fitter.data.loaded_spec_data["spectrum1"]["count_rate"])
axes[0].set_ylim([1e-1, y_max_lim])
axes[0].set_xlim([2, 7])
plt.tight_layout()
plt.show()
