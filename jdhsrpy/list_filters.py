import astropy.units as u
import numpy as np
from astropy.time import Time

from jdhsrpy import IMAGE_RES, NUSTAR_EPOCH

__all__ = [
    "by_detector",
    "by_energy",
    "by_good_pix",
    "by_gradezero",
    "by_region",
    "by_time",
    "event_filter",
]


def by_good_pix(event_data, fpm):
    """Do some basic filtering on known bad pixels.

    Parameters
    ----------
    event_data: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.

    fpm: {"A" | "B"}
        Which FPM you're filtering on. Assumes A if not set.

    Returns
    -------

    goodinds: iterable
        Index of event_data that passes the filtering.
    """
    # Hot pixel filters

    # FPMA or FPMB

    if fpm.find("B") == -1:
        pix_filter = np.invert(
            (event_data["DET_ID"] == 2)
            & (event_data["RAWX"] == 16)
            & (event_data["RAWY"] == 5)
            | (event_data["DET_ID"] == 2)
            & (event_data["RAWX"] == 24)
            & (event_data["RAWY"] == 22)
            | (event_data["DET_ID"] == 2)
            & (event_data["RAWX"] == 27)
            & (event_data["RAWY"] == 6)
            | (event_data["DET_ID"] == 2)
            & (event_data["RAWX"] == 27)
            & (event_data["RAWY"] == 21)
            | (event_data["DET_ID"] == 3)
            & (event_data["RAWX"] == 22)
            & (event_data["RAWY"] == 1)
            | (event_data["DET_ID"] == 3)
            & (event_data["RAWX"] == 15)
            & (event_data["RAWY"] == 3)
            | (event_data["DET_ID"] == 3)
            & (event_data["RAWX"] == 5)
            & (event_data["RAWY"] == 5)
            | (event_data["DET_ID"] == 3)
            & (event_data["RAWX"] == 22)
            & (event_data["RAWY"] == 7)
            | (event_data["DET_ID"] == 3)
            & (event_data["RAWX"] == 16)
            & (event_data["RAWY"] == 11)
            | (event_data["DET_ID"] == 3)
            & (event_data["RAWX"] == 18)
            & (event_data["RAWY"] == 3)
            | (event_data["DET_ID"] == 3)
            & (event_data["RAWX"] == 24)
            & (event_data["RAWY"] == 4)
            | (event_data["DET_ID"] == 3)
            & (event_data["RAWX"] == 25)
            & (event_data["RAWY"] == 5)
        )
    else:
        pix_filter = np.invert(
            (event_data["DET_ID"] == 0)
            & (event_data["RAWX"] == 24)
            & (event_data["RAWY"] == 24)
        )

    inds = (pix_filter).nonzero()
    goodinds = inds[0]

    return event_data[goodinds]


def by_energy(event_data, energy_low=2.5, energy_high=10):
    """Apply energy filtering to the data.

    Parameters
    ----------
    event_data: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.

    energy_low: float
        Low-side energy bound for the map you want to produce (in keV).
        Defaults to 2.5 keV.

    energy_high: float
        High-side energy bound for the map you want to produce (in keV).
        Defaults to 10 keV.
    """
    pilow = (energy_low - 1.6) / 0.04
    pihigh = (energy_high - 1.6) / 0.04
    pi_filter = (event_data["PI"] >= pilow) & (event_data["PI"] < pihigh)
    inds = (pi_filter).nonzero()
    goodinds = inds[0]

    return event_data[goodinds]


def by_gradezero(event_data):
    """Only accept counts with GRADE==0.

    Parameters
    ----------
    event_data: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.

    Returns
    -------
    goodinds: iterable
        Index of event_data that passes the filtering.
    """
    # Grade filter
    return by_grade(event_data, 0)


def by_grade(event_data, grade: int):
    """Only accept counts with GRADE==`grade`.

    Parameters
    ----------
    event_data: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.

    grade : `int`
        Must be a valid NuSTAR event grade.

    Returns
    -------
    goodinds: iterable
        Index of event_data that passes the filtering.
    """

    # Grade filter
    grade_filter = event_data["GRADE"] == grade
    inds = (grade_filter).nonzero()
    goodinds = inds[0]

    return event_data[goodinds]


def by_time(event_data, start, end):
    """Only include counts within a given time range.

    Parameters
    ----------
    event_data: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.

    start : `str`
        Input two times in the form 'yyyy/mm/dd, HH:MM:SS'
        (e.g. '2019/03/06, 16:45:30') for the time range,
        default is the whole observation for the file.

    Returns
    -------
    goodinds: iterable
        Index of event_data that lies within the time range.
    """

    tstart = Time(start, format="isot", scale="utc")
    tend = Time(end, format="isot", scale="utc")
    tstart_s = (
        tstart - NUSTAR_EPOCH
    ).sec  # both dates are converted to number of seconds from 2010-Jan-1
    tend_s = (tend - NUSTAR_EPOCH).sec
    tmrng = [tstart_s, tend_s]

    time_filter = (event_data["TIME"] > tmrng[0]) & (event_data["TIME"] < tmrng[1])
    inds = (time_filter).nonzero()
    goodinds = inds[0]

    return event_data[goodinds]


def event_filter(event_data, fpm="FPMA", energy_low=2.5, energy_high=10, tmrng=None):
    # was event_filter(event_data, fpm='FPMA', energy_low=2.5, energy_high=10) # Kris #
    """All in one filter module. By default applies an energy cut,
        selects only events with grade == 0, and removes known hot pixel.

        Note that this module returns a cleaned eventlist rather than
        the indices to the cleaned events.

    Parameters
    ----------
    event_data: FITS data structure
        This should be an hdu.data structure from a NuSTAR FITS file.

    fpm: {"FPMA" | "FPMB"}
        Which FPM you're filtering on. Defaults to FPMA.

    energy_low: float
        Low-side energy bound for the map you want to produce (in keV).
        Defaults to 2.5 keV.

    energy_high: float
        High-side energy bound for the map you want to produce (in keV).
        Defaults to 10 keV.

    tmrng : list of strings, length 2
        Input two times in the form 'yyyy-mm-ddTHH:MM:SS.ZZZ'
        (e.g. '2010-01-01T00:00:00.000') for the time range,
        default is the whole observation for the file.

    Returns
    -------

    cleanevt: FITS data class.
        This is the subset of event_data that pass the data selection cuts.
    """

    evt_timefilter = (
        by_time(event_data, tmrng[0], tmrng[1]) if tmrng is not None else event_data
    )
    evt_badfilter = by_good_pix(evt_timefilter, fpm=fpm)
    evt_energy = by_energy(
        evt_badfilter, energy_low=energy_low, energy_high=energy_high
    )
    return by_gradezero(evt_energy)


def by_detector(event_data, detector):
    """Filter and return event list filtered for the desired detector."""
    if detector not in range(4):
        raise ValueError("the `detector` input much be an int in [0,1,2,3].")
    return event_data[event_data["DET_ID"] == detector]


def by_region(event_data, bottom_left, top_right):
    """Filter and return event list filtered for a square region."""
    pix_map_centre = 1500
    sol_x_pix = (event_data["X"] - pix_map_centre) << u.pixel
    sol_y_pix = (event_data["Y"] - pix_map_centre) << u.pixel
    sol_x_arc = sol_x_pix * IMAGE_RES
    sol_y_arc = sol_y_pix * IMAGE_RES
    reg_filter = (
        (bottom_left[0] <= sol_x_arc)
        & (sol_x_arc < top_right[0])
        & (bottom_left[1] <= sol_y_arc)
        & (sol_y_arc < top_right[1])
    )
    inds = (reg_filter).nonzero()[0]
    return event_data[inds]
