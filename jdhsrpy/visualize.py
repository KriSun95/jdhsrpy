from datetime import datetime

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import time_support
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.polygon import Polygon

__all__ = [
    "chu_plot",
    "livetime_plot",
    "optical_axis_plot",
    "spectrum_plot",
    "time_profile_plot",
    "vertical_line_of_time",
]


def time_profile_plot(times, ydata, axes=None, **kwargs):
    """Plot a time profile (on a given axes if given).

    Parameters
    ----------
    times : `~astropy.time.Time`|`~datetime.datetime`
        The time edges.

    ydata : `list`|`~numpy.ndarray`
        The data varying with time. Length will be one less the `times`.
        length.

    axes : axis object
        The axis to be plotted on. If None then ``plt.gca()`` is used.
        Default: None

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    time_support(format="unix_tai")
    axes.stairs(ydata, to_plotting_datetimes(times), **kwargs)
    return axes


def to_plotting_datetimes(times):
    """Convert times to ones that plot well."""
    if isinstance(times, Time):
        return Time(times, format="unix_tai", scale="utc").datetime
    elif isinstance(times, datetime):
        return times
    else:
        raise TypeError("The `times` variable is not supported.")


def vertical_line_of_time(time: str | Time | datetime, axes=None, **kwargs):
    """Use `axvline` to plot `time` input.

    Parameters
    ----------
    time : `str`|`astropy.time.core.Time`|`datetime`
        If `str` then need ISOT format. E.g., 2021-11-20T02:25:30.

    axes : axis object
        The axis to be plotted on. If None then ``plt.gca()`` is used.
        Default: None

    Returns
    -------
    Axes input.
    """
    axes = plt if axes is None else axes
    if isinstance(time, str):
        converted_time = Time(time, format="isot", scale="utc").datetime
    else:
        converted_time = to_plotting_datetimes(time)

    axes.axvline(x=converted_time, **kwargs)
    return axes


def chu_plot(times, ydata, label_map, axes=None, **kwargs):
    """Plot a time profile (on a given axes if given).

    Parameters
    ----------
    times : `~astropy.time.Time`|`~datetime.datetime`
        The time centers.

    ydata : `list`|`~numpy.ndarray`
        The data varying with time. Length will be the same as `times`.
        length.

    axes : axis object
        The axis to be plotted on. If None then ``plt.gca()`` is used.
        Default: None

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    time_support(format="unix_tai")
    axes.plot(to_plotting_datetimes(times), ydata, "x", **kwargs)
    axes.axes.set_yticks(list(label_map.values()))
    axes.axes.set_yticklabels(list(label_map.keys()))
    return axes


def livetime_plot(times, ydata, axes=None, **kwargs):
    """Plot a time profile (on a given axes if given).

    Parameters
    ----------
    times : `~astropy.time.Time`|`~datetime.datetime`
        The time centers.

    ydata : `list`|`~numpy.ndarray`
        The data varying with time. Length will be the same as `times`.
        length.

    axes : axis object
        The axis to be plotted on. If None then ``plt.gca()`` is used.
        Default: None

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    _defaults = {
        "drawstyle": "steps-mid",
    }
    _plot_settings = _defaults | kwargs
    plt.semilogy(to_plotting_datetimes(times), ydata, **_plot_settings)
    return axes


def spectrum_plot(energy_bins, ydata, ydata_err, axes=None, **kwargs):
    """Plot a spectrum profile (on a given axes if given).

    Parameters
    ----------
    energy_bins : `~astropy.unit.Quantity`
        The energy bin edges. Should be one more than the `ydata` length.

    ydata :  `~astropy.unit.Quantity`
        The spectrum.

    axes : axis object
        The axis to be plotted on. If None then ``plt.gca()`` is used.
        Default: None

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    _defaults = {
        "color": "k",
        "fmt": ".",
        "markersize": 0.01,
    }
    _plot_settings = _defaults | kwargs

    energy_mids = (energy_bins[:-1] + energy_bins[1:]) / 2
    energy_half_width = (energy_bins[1:] - energy_bins[:-1]) / 2

    ## plot spectrum created from Python
    axes.errorbar(
        energy_mids, ydata, xerr=energy_half_width, yerr=ydata_err, **_plot_settings
    )

    axes.set_ylabel(f"{ydata.unit:latex}")
    axes.tick_params(axis="both", which="both", direction="in")
    axes.set_yscale("log")
    axes.set_xlabel(f"Energy [{energy_bins.unit:latex}]")
    return axes


def optical_axis_plot(times, positions, axes=None, **kwargs):
    """Plot a time profile (on a given axes if given).

    Parameters
    ----------
    times : `~astropy.time.Time`|`~datetime.datetime`
        The time centers.

    positions : `list`|`~numpy.ndarray`|`astropy.units.Quantity`
        The optical axis positions for either the X or Y axis.

    axes : axis object
        The axis to be plotted on. If None then ``plt.gca()`` is used.
        Default: None

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    _defaults = {
        "drawstyle": "steps-mid",
    }
    _plot_settings = _defaults | kwargs
    plt.plot(to_plotting_datetimes(times), positions, **_plot_settings)
    return axes


def shape_plot(
    shapely_polygon, coordinate_frame, axes=None, sky_coord_kwargs=None, **kwargs
):
    """Plot a ``shapely.geometry.polygon.Polygon``.

    Parameters
    ----------
    shapely_polygon :  `shapely.geometry.polygon.Polygon` | `geometry.multipolygon.MultiPolygon`
        Will assume ``shapely_polygon.exterior.coords`` is in arcseconds
        by default (see `sky_coord_kwargs` for more). If ``geometry.multipolygon.MultiPolygon``
        (say a symmetric different has been used, i.e. multiple separate
        polygons exist) then these will be plotted with colours cycling
        through Matplotib's color cycle. Although, passing a line
        property such as ``c`` or ``color`` as a ``kwargs`` will apply
        to all polygons.

    coordinate_frame : Sunpy coordinate map frame

    axes : axis object
        The axis to be plotted on. If None then ``plt.gca()`` is used.
        Default: None

    sky_coord_kwargs : `dict` | `None
        Passed to ``astropy.coordinates.SkyCoord`` If None then the
        default unit for the ``shapely.geometry.polygon.Polygon`` will
        be arcseconds.
        Default: None

    **kwargs :
        Passed to ``plot_coord``.

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    sky_coord_kwargs = {} if sky_coord_kwargs is None else sky_coord_kwargs
    _sky_coord_kwargs = {"unit": u.arcsec} | sky_coord_kwargs

    if isinstance(shapely_polygon, Polygon):
        vertices = np.array(shapely_polygon.exterior.coords)
        points_plot(
            vertices[:, 0],
            vertices[:, 1],
            coordinate_frame,
            axes=axes,
            sky_coord_kwargs=_sky_coord_kwargs,
            **kwargs,
        )

    elif isinstance(shapely_polygon, MultiPolygon):
        for shape in shapely_polygon.geoms:
            vertices = np.array(shape.exterior.coords)
            points_plot(
                vertices[:, 0],
                vertices[:, 1],
                coordinate_frame,
                axes=axes,
                sky_coord_kwargs=_sky_coord_kwargs,
                **kwargs,
            )
    else:
        raise TypeError(
            "Input for ``shapely_polygon`` is not of type ``geometry.multipolygon.MultiPolygon`` or ``geometry.polygon.Polygon``. "
        )
    return axes


def points_plot(x, y, coordinate_frame, axes=None, sky_coord_kwargs=None, **kwargs):
    """Plot a set of coordinates.

    Parameters
    ----------
    x, y :  `list` | array-like
        The x and y coordinates to be plotted. Given in solar pixel
        coordinates and plotted in helioprojective.

    coordinate_frame : Sunpy coordinate map frame

    axes : axis object
        The axis to be plotted on. If None then ``plt.gca()`` is used.
        Default: None

    sky_coord_kwargs : `dict` | `None
        Passed to ``astropy.coordinates.SkyCoord`` If None then the
        default unit for the ``shapely.geometry.polygon.Polygon`` will
        be arcseconds.
        Default: None

    **kwargs :
        Passed to ``plot_coord``.

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    sky_coord_kwargs = {} if sky_coord_kwargs is None else sky_coord_kwargs
    _sky_coord_kwargs = {"unit": u.arcsec} | sky_coord_kwargs
    sc = SkyCoord(x, y, frame=coordinate_frame, **_sky_coord_kwargs)
    axes.plot_coord(sc, **kwargs)
    return axes
