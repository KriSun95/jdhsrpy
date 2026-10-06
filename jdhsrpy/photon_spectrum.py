import warnings

import astropy.units as u
import numpy as np
from numpy.exceptions import VisibleDeprecationWarning
from sunkit_spex.legacy.emission import bremsstrahlung_thick_target
from sunkit_spex.legacy.thermal import thermal_emission

warnings.filterwarnings("ignore", category=RuntimeWarning)
try:
    warnings.filterwarnings("ignore", category=VisibleDeprecationWarning)
except AttributeError:
    warnings.filterwarnings("ignore", category=VisibleDeprecationWarning)

__all__ = ["thermal", "thick_target"]


def thick_target(energies, total_eflux, index, e_c):
    """Single electron power-law thick-target emission

    Parameters
    ---------
    energies : `~numpy.ndarray`
        A 1D array of the energy bin edges.

    total_eflux : `int` | `float`
        The electron flux in 1e35 e/s.

    index : `int` | `float`
        The power-law index.

    e_c : `int` | `float`
        The low-energy cut-off in keV.

    Returns
    -------
    : `~numpy.ndarray`
        The photon spectrum in ph/keV/s/cm^2 with length one less than
        the `energy` array.
    """
    energies = energies.value
    high_break = np.float64(energies.max() * 10)
    output = bremsstrahlung_thick_target(
        photon_energies=energies,
        p=index,
        eebrk=high_break,
        q=20,
        eelow=e_c,
        eehigh=high_break,
    )

    output[np.isnan(output)] = 0
    output[~np.isfinite(output)] = 0

    # convert to 1e35 e-/s
    return (output * total_eflux * 1e35) << u.ph / u.keV / u.s / u.cm**2


def thermal(energies, temperature, emission_measure):
    """Thermal bremsstrahlung emission

    Parameters
    ---------
    energies : `~numpy.ndarray`
        A 1D array of the energy bin edges.

    temperature : `int` | `float`
        The plasma temperature in MK.

    emission_measure : `int` | `float`
        The plasma emission measure in cm^-3.

    Returns
    -------
    : `~numpy.ndarray`
        The photon spectrum in ph/keV/s/cm^2 with length one less than
        the `energy` array.
    """
    return thermal_emission(energies, temperature, emission_measure)
