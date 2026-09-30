
import astropy
import astropy.units as u
import numpy as np

from jdhsr.io import read_nustar_pha 

__all__ = ["get_observable_info",
           "get_effective_area_info", 
           "get_response_info",
           "nustar_pha_spectrum",
           ]

def get_observable_info(pha_data:astropy.io.fits.fitsrec.FITS_rec, pha_header:astropy.io.fits.fitsrec.FITS_rec):
    """Extract the channel, observable, and livetime from NuSTAR PHA file."""
    return pha_data["channel"]<<u.dimensionless_unscaled, pha_data["counts"]<<u.ct, pha_header["LIVETIME"]<<u.second

def get_effective_area_info(arf_data:astropy.io.fits.fitsrec.FITS_rec):
    """Extract the channel, observable, and livetime from NuSTAR ARF file."""
    return arf_data["energ_lo"]<<u.keV, arf_data["energ_hi"]<<u.keV, arf_data["specresp"]<<u.cm**2

def get_response_info(rmf_cdata:astropy.io.fits.fitsrec.FITS_rec, rmf_pdata:astropy.io.fits.fitsrec.FITS_rec):
    """Extract the channel, observable, and livetime from NuSTAR RMF file."""
    return (rmf_cdata["channel"]<<u.dimensionless_unscaled, rmf_cdata["e_min"]<<u.keV, rmf_cdata["e_max"]<<u.keV), (rmf_pdata["energ_lo"]<<u.keV, rmf_pdata["energ_hi"]<<u.keV, rmf_pdata["n_grp"]<<u.dimensionless_unscaled, rmf_pdata["f_chan"], rmf_pdata["n_chan"], rmf_pdata["matrix"])

def standard_spectrum_axis():
    """Defines the usual, native NuSTAR energy binning."""
    _standard_energy_start = 1.6 << u.keV
    _standard_num_of_channels = 4096
    _standard_energy_binning = 0.04 << u.keV
    _standard_max_energy = _standard_num_of_channels*_standard_energy_binning + _standard_energy_start
    e_lo = np.arange(_standard_energy_start.value, _standard_max_energy.value, _standard_energy_binning.value) << _standard_energy_start.unit
    e_end = e_lo[-1] + _standard_energy_binning
    return np.append(e_lo, e_end)

def nustar_pha_spectrum(file, counts_only=False):
    """Takes a .pha file and returns plotting information.
    
    Parameters
    ----------
    file : Str
        String for the .pha file of the spectrum under investigation.

    counts_only : `bool`
        If True, then return just the counts spectrum, else return the 
        spectrum in counts/keV/s.
        Default: False
            
    Returns
    -------
    The energy bin edges, count or count rate per keV, and its error. 
    """

    _, counts, livetime = get_observable_info(*read_nustar_pha(file))
    
    e_bins = standard_spectrum_axis()
    e_bin_size = np.diff(e_bins)[0]

    if counts_only:
        livetime = 1
        e_bin_size = 1
    
    cts = (counts / e_bin_size) / livetime # now in cts keV^-1 s^-1
    cts_err = (np.sqrt(counts.value)<<counts.unit / e_bin_size) / livetime
    
    return e_bins, cts, cts_err