import numpy as np
from astropy.time import Time

from jdhsrpy import NUSTAR_EPOCH

__all__ = [
    "by_detector",
    "by_energy",
    "by_good_pix",
    "by_gradezero",
    "by_time",
    "event_filter",
]


def by_good_pix(evtdata, fpm):
    """Do some basic filtering on known bad pixels.

    Parameters
    ----------
    evtdata: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.

    fpm: {"A" | "B"}
        Which FPM you're filtering on. Assumes A if not set.

    Returns
    -------

    goodinds: iterable
        Index of evtdata that passes the filtering.
    """
    # Hot pixel filters

    # FPMA or FPMB

    if fpm.find("B") == -1:
        pix_filter = np.invert(
            (evtdata["DET_ID"] == 2) & (evtdata["RAWX"] == 16) & (evtdata["RAWY"] == 5)
            | (evtdata["DET_ID"] == 2)
            & (evtdata["RAWX"] == 24)
            & (evtdata["RAWY"] == 22)
            | (evtdata["DET_ID"] == 2)
            & (evtdata["RAWX"] == 27)
            & (evtdata["RAWY"] == 6)
            | (evtdata["DET_ID"] == 2)
            & (evtdata["RAWX"] == 27)
            & (evtdata["RAWY"] == 21)
            | (evtdata["DET_ID"] == 3)
            & (evtdata["RAWX"] == 22)
            & (evtdata["RAWY"] == 1)
            | (evtdata["DET_ID"] == 3)
            & (evtdata["RAWX"] == 15)
            & (evtdata["RAWY"] == 3)
            | (evtdata["DET_ID"] == 3) & (evtdata["RAWX"] == 5) & (evtdata["RAWY"] == 5)
            | (evtdata["DET_ID"] == 3)
            & (evtdata["RAWX"] == 22)
            & (evtdata["RAWY"] == 7)
            | (evtdata["DET_ID"] == 3)
            & (evtdata["RAWX"] == 16)
            & (evtdata["RAWY"] == 11)
            | (evtdata["DET_ID"] == 3)
            & (evtdata["RAWX"] == 18)
            & (evtdata["RAWY"] == 3)
            | (evtdata["DET_ID"] == 3)
            & (evtdata["RAWX"] == 24)
            & (evtdata["RAWY"] == 4)
            | (evtdata["DET_ID"] == 3)
            & (evtdata["RAWX"] == 25)
            & (evtdata["RAWY"] == 5)
        )
    else:
        pix_filter = np.invert(
            (evtdata["DET_ID"] == 0) & (evtdata["RAWX"] == 24) & (evtdata["RAWY"] == 24)
        )

    inds = (pix_filter).nonzero()
    goodinds = inds[0]

    return evtdata[goodinds]


def by_energy(evtdata, energy_low=2.5, energy_high=10):
    """Apply energy filtering to the data.

    Parameters
    ----------
    evtdata: FITS data class
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
    pi_filter = (evtdata["PI"] >= pilow) & (evtdata["PI"] < pihigh)
    inds = (pi_filter).nonzero()
    goodinds = inds[0]

    return evtdata[goodinds]


def by_gradezero(evtdata):
    """Only accept counts with GRADE==0.

    Parameters
    ----------
    evtdata: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.

    Returns
    -------
    goodinds: iterable
        Index of evtdata that passes the filtering.
    """
    # Grade filter
    return by_grade(evtdata, 0)


def by_grade(evtdata, grade: int):
    """Only accept counts with GRADE==`grade`.

    Parameters
    ----------
    evtdata: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.

    grade : `int`
        Must be a valid NuSTAR event grade.

    Returns
    -------
    goodinds: iterable
        Index of evtdata that passes the filtering.
    """

    # Grade filter
    grade_filter = evtdata["GRADE"] == grade
    inds = (grade_filter).nonzero()
    goodinds = inds[0]

    return evtdata[goodinds]


def by_time(evtdata, tmrng):
    """Only include counts within a given time range.

    Parameters
    ----------
    evtdata: FITS data class
        This should be an hdu.data structure from a NuSTAR FITS file.

    tmrng : list of length 2
        Input two times in the form 'yyyy/mm/dd, HH:MM:SS'
        (e.g. '2019/03/06, 16:45:30') for the time range,
        default is the whole observation for the file.

    Returns
    -------
    goodinds: iterable
        Index of evtdata that lies within the time range.
    """
    if tmrng is None:
        return np.arange(len(evtdata))

    tstart = Time(tmrng[0], format="isot", scale="utc")
    tend = Time(tmrng[1], format="isot", scale="utc")
    tstart_s = (
        tstart - NUSTAR_EPOCH
    ).sec  # both dates are converted to number of seconds from 2010-Jan-1
    tend_s = (tend - NUSTAR_EPOCH).sec
    tmrng = [tstart_s, tend_s]

    time_filter = (evtdata["TIME"] > tmrng[0]) & (evtdata["TIME"] < tmrng[1])
    inds = (time_filter).nonzero()
    goodinds = inds[0]

    return evtdata[goodinds]


def event_filter(evtdata, fpm="FPMA", energy_low=2.5, energy_high=10, tmrng=None):
    # was event_filter(evtdata, fpm='FPMA', energy_low=2.5, energy_high=10) # Kris #
    """All in one filter module. By default applies an energy cut,
        selects only events with grade == 0, and removes known hot pixel.

        Note that this module returns a cleaned eventlist rather than
        the indices to the cleaned events.

    Parameters
    ----------
    evtdata: FITS data structure
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
        This is the subset of evtdata that pass the data selection cuts.
    """

    evt_timefilter = by_time(evtdata, tmrng)
    evt_badfilter = by_good_pix(evt_timefilter, fpm=fpm)
    evt_energy = by_energy(
        evt_badfilter, energy_low=energy_low, energy_high=energy_high
    )
    return by_gradezero(evt_energy)


def by_detector(evtdata, detector):
    """Filter and return event list filtered for the desired detector."""
    if detector not in range(4):
        raise ValueError("the `detector` input much be an int in [0,1,2,3].")
    return evtdata[evtdata["DET_ID"] == detector]
