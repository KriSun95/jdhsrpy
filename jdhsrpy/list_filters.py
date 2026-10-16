import astropy.units as u
import numpy as np
from astropy.time import Time
from shapely import points

from jdhsrpy import NUSTAR_EPOCH
from jdhsrpy.shapes import circle_shape, rectangle_shape
from jdhsrpy.utils import pixel_coords_to_solar_coords

__all__ = [
    "by_basic_rectangle_region",
    "by_circle_region",
    "by_detector",
    "by_energy",
    "by_general_region",
    "by_good_pix",
    "by_grade",
    "by_gradezero",
    "by_rectangle_region",
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


def by_basic_rectangle_region(event_data, bottom_left, top_right):
    """Filter and return event list filtered for a rectangle region.

    This existed before the use of the ``shapely`` module where the
    region selection just simply checks the count location against the
    X and Y bounds of the box.

    This function is kept for legacy reasons and can serve as a nice
    check for the other region event list filters.

    Parameters
    ----------
    event_data: FITS data structure
        This should be an hdu.data structure from a NuSTAR FITS file.

    bottom_left : `~astropy.units.Quantity`
        The bottom left corner (X, Y) of the desired box. Should be to
        arcseconds.

    top_right : `~astropy.units.Quantity`
        The top right corner (X, Y) of the desired box. Should be to
        arcseconds.
    """
    sol_x_arc = pixel_coords_to_solar_coords(event_data["X"])
    sol_y_arc = pixel_coords_to_solar_coords(event_data["Y"])
    reg_filter = (
        (bottom_left[0] <= sol_x_arc)
        & (sol_x_arc < top_right[0])
        & (bottom_left[1] <= sol_y_arc)
        & (sol_y_arc < top_right[1])
    )
    inds = (reg_filter).nonzero()[0]
    return event_data[inds]


def by_general_region(event_data, shape):
    """Filter and return event list filtered for a general region.

    Parameters
    ----------
    event_data: FITS data structure
        This should be an hdu.data structure from a NuSTAR FITS file.

    shape : `~shapely.geometry.polygon.Polygon`
    """
    # enforce everything in arcseconds before operations
    sol_x_arc = pixel_coords_to_solar_coords(event_data["X"]) << u.arcsec
    sol_y_arc = pixel_coords_to_solar_coords(event_data["Y"]) << u.arcsec

    map_points = points(np.vstack((sol_x_arc.value, sol_y_arc.value)).T)
    in_shape = np.array([(mp.within(shape) or mp.touches(shape)) for mp in map_points])
    return event_data[in_shape]


def by_rectangle_region(event_data, bottom_left, top_right):
    """Filter and return event list filtered for a rectangle region.

    Is identical to the ``by_basic_rectangle_region`` function but uses
    the far more powerful ``shapely`` Python package.

    Parameters
    ----------
    event_data: FITS data structure
        This should be an hdu.data structure from a NuSTAR FITS file.

    bottom_left : `~astropy.units.Quantity`
        The bottom left corner (X, Y) of the desired box. Should be to
        arcseconds.

    top_right : `~astropy.units.Quantity`
        The top right corner (X, Y) of the desired box. Should be to
        arcseconds.
    """
    shape = rectangle_shape(bottom_left, top_right)
    return by_general_region(event_data, shape)


def by_circle_region(event_data, center, radius):
    """Filter and return event list filtered for a circle region.

    Parameters
    ----------
    event_data: FITS data structure
        This should be an hdu.data structure from a NuSTAR FITS file.

    center : `~astropy.units.Quantity`
        Thecenter (X, Y) of the desired circle. Should be to arcseconds.

    radius : `~astropy.units.Quantity`
        The radius R of the desired circle. Should be to arcseconds.
    """
    shape = circle_shape(center, radius)
    return by_general_region(event_data, shape)
