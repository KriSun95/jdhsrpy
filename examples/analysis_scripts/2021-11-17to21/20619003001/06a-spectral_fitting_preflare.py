import ntpath
import os

import matplotlib.pyplot as plt
import numpy as np

from jdhsrpy import spectral_fitting

DIRECTORY, FILENAME = ntpath.split(__file__)
DESCRIPTION = f"ISOTHERMAL_PREFLARE_FROM_{FILENAME[:-3]}"

SAVE_DIRECTORY = os.path.join(DIRECTORY, DESCRIPTION)
os.makedirs(SAVE_DIRECTORY, exist_ok=True)
SAVE_PLOTS = False

# define the spectral files (ARF and RMF will be found automatically)
nustar_pha_files = {"FPMA":"/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/event_cl/20211120022210_to_20211120022530/nu20619003001A06_cl_grade0_sr.pha",
                    "FPMB":"/Users/kris/Documents/umnPostdoc/projects/analysis/nustarNov2021/data/nsNov2021on17-19-21/nustarFiles/nsNov19/20619003001/event_cl/20211120022210_to_20211120022530/nu20619003001B06_cl_grade0_sr.pha",
                    }
# pass the files to the fitter object
fitter = spectral_fitting.get_fitter_object(*nustar_pha_files.values())
# define the model we want
spectral_fitting.set_single_thermal_model(fitter)
# Define the energy range for the fit in keV
fitter.energy_fitting_range = [2.5, 6.5]
# Give starting points and bounds for the parameters
fitter.params["T1_spectrum1"] = {"Value":3, "Bounds":(1.1, 5)}
fitter.params["EM1_spectrum1"] = {"Value":5.5e-2, "Bounds":(5e-2, 5e0)}
fitter.params["C_spectrum1"] = "freeze"
fitter.params["C_spectrum2"] = {"Status":"free", "Value":1, "Bounds":(0.8, 1.1)}

# fit by optimisation
fitter.fit()
# plot the result
plt.figure(figsize=(25,7))
axes, res_axes = fitter.plot()
y_max_lim = 1.2*np.max(fitter.data.loaded_spec_data["spectrum1"]["count_rate"])
for axis in axes:
    axis.set_ylim([1e-1, y_max_lim])
    axis.set_xlim([2, 12])
dets = list(nustar_pha_files.keys())
axes[0].set_title(dets[0])
axes[1].set_title(dets[1])
plt.tight_layout()
if SAVE_PLOTS:
    plt.savefig(os.path.join(SAVE_DIRECTORY, f"{DESCRIPTION}_minimize_fit.png"))
plt.show()

# fit via MCMC
mcmc_result = fitter.run_mcmc(steps_per_walker=1_000)
fitter.burn_mcmc = 500
# log-probability chain as the walkers wander
plt.figure()
fitter.plot_log_prob_chain()
plt.tight_layout()

if SAVE_PLOTS:
    plt.savefig(os.path.join(SAVE_DIRECTORY, f"{DESCRIPTION}_MCMC_prob_chain.png"))
plt.show()
# corner plot for the parameters
corner_plot = fitter.corner_mcmc()
plt.tight_layout()

if SAVE_PLOTS:
    plt.savefig(os.path.join(SAVE_DIRECTORY, f"{DESCRIPTION}_MCMC_corner.png"))
plt.show()
# plot the resulting fit
plt.figure(figsize=(25,7))
axes, res_axes = fitter.plot()
for axis in axes:
    axis.set_ylim([1e-1, y_max_lim])
    axis.set_xlim([2, 12])
dets = list(nustar_pha_files.keys())
axes[0].set_title(dets[0])
axes[1].set_title(dets[1])
plt.tight_layout()

if SAVE_PLOTS:
    plt.savefig(os.path.join(SAVE_DIRECTORY, f"{DESCRIPTION}_MCMC_fit.png"))
plt.show()
