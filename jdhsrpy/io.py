"""
Module containing functions that load HEASARC compliant files.
"""

import warnings

import astropy.units as u
from astropy.io import fits

from jdhsrpy.utils import (
    create_region_file_text,
    deconstruct_degminsec,
    deconstruct_hrminsec,
)

__all__ = [
    "read_heasarc_arf",
    "read_heasarc_rmf",
    "read_nustar_pha",
    "read_region_file",
    "write_region_file",
]


def read_nustar_pha(file: str):
    """
    Read a `.pha` file and extract data ande header information.

    Parameters
    ----------
    file : `str`, `file-like` or `pathlib.Path`
        A `.pha` file (see `~astropy.fits.io.open` for details).

    Returns
    -------
    The counts and channel data component and the livetime header
    component.
    """
    with fits.open(file) as hdul:
        data = hdul[1].data
        header = hdul[0].header

    return data, header


def read_heasarc_arf(file: str):
    """
    Read a HEASARC compliant `.arf` file and extract useful information.

    [1] https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html

    Parameters
    ----------
    file :  `str`, `file-like` or `pathlib.Path`
        A `.arf` file (see `~astropy.fits.io.open` for details ).

    Returns
    -------
    The effective area data component from the file.
    """
    data = None
    with fits.open(file) as hdul:
        for hdu in hdul:
            hdu_contents = hdu.header.get("HDUCLAS2", None)
            if hdu_contents == "SPECRESP":
                data = hdu.data

    return data


def read_heasarc_rmf(file: str):
    """
    Read a HEASARC compliant `.rmf` file and extract useful information.

    [1] https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html

    Parameters
    ----------
    file :  `str`, `file-like` or `pathlib.Path`
        A `.rmf` file (see `~astropy.fits.io.open` for details).

    Returns
    -------
    The channel bin data mapping (channel number to energy for the
    corresponding PHA file) and the RMF & photon channel data.
    """
    channel_data = None
    rmf_and_photon_data = None
    with fits.open(file) as hdul:
        for hdu in hdul:
            hdu_contents = hdu.header.get("HDUCLAS2", None)
            if hdu_contents == "EBOUNDS":
                channel_data = hdu.data
            elif hdu_contents == "RSP_MATRIX":
                rmf_and_photon_data = hdu.data

    return channel_data, rmf_and_photon_data


def read_filter_file(file: str):
    """Read the NuSTAR filter file (ext. MKF) and return the roool angle.

    Times for each roll angle are also returned.

    Parameters
    ----------
    file :  `str`, `file-like` or `pathlib.Path`
        A `.mkf` file.
    """
    with fits.open(file) as hdul:
        filter_data = hdul[1].data
    return filter_data["TIME"] << u.s, filter_data["ROLL"] << u.deg


def read_region_file(reg_file):
    """Get the RA, Dec, and radius of the circle region in the REG file.

    The REG file should have been saved using FK5 or ICRS as its WCS.

    The circle center (ra, dec) can be in degrees or sexagesimal.

    Parameters
    ----------
    reg_file : str
        A region file to be read. For NuSTAR, this should be a file with
        something like \'circle(ra, dec, rad)\' in it somewhere.

    Returns
    -------
    The RA, Dec, and redius of the circle region.
    """

    with open(reg_file) as f:
        contents = f.read()

    if ("fk5" not in contents) and ("icrs" not in contents):
        warnings.warn("The region file is not in FK5 or ICRS coordinates.")
        return

    _ra, _dec, _rad_arcsec = contents.split("(")[-1].split(")")[0].split(",")
    if ":" in _ra:
        ra = deconstruct_hrminsec(_ra)
        dec = deconstruct_degminsec(_dec)
    else:
        ra = float(_ra) << u.deg
        dec = float(_dec) << u.deg
    rad_arcsec = float(_rad_arcsec[:-1]) << u.arcsec

    return ra, dec, rad_arcsec


def write_region_file(ra, dec, radius, filename, frame="fk5"):
    """Create the text to be saved in th region file.

    The body of the text has been taken from an existing REG file from
    SAOImageDS9 and then we replace the relevant sections with the input
    values.

    Parameters
    ----------
    ra : `astropy.units.Quantity`
        The right ascension of circle centre in, or convertible to,
        degrees.

    dec : `astropy.units.Quantity`
        The declination of circle centre in, or convertible to, degrees.

    radius : `astropy.units.Quantity`
        The radius of circle in, or convertible to, degrees.

    filename : `str`
        The name of the output region file. If it does not end in
        \".reg\" then this will be appended.

    frame : `str`
        The coordinate frame. In the NuSTAR file header, it is stated
        that the frame is FK5 so that is the default. However,
        SAOImageDS9 can also output the frame as ICRS (I have not seen
        the difference in RA, Dec, or radius values between them yet).
        Therefore, the options here are \"fk5\" and \"icrs\" (note the
        lower case).
    """
    filename = filename if filename.endswith(".reg") else f"{filename}.reg"
    contents = create_region_file_text(ra, dec, radius, frame=frame)

    with open(filename, "w") as f:
        f.write(contents)
