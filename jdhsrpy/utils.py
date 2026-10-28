import re

import astropy.units as u
import numpy as np
from shapely.geometry.polygon import Polygon
from sunpy.coordinates import sun

from jdhsrpy import IMAGE_COLORMAP, IMAGE_RES, IMAGE_SCALE

__all__ = [
    "arcsec_points_to_shape",
    "assign_plot_settings",
    "deconstruct_degminsec",
    "only_numbers",
    "pixel_coords_to_solar_coords",
    "radec_to_helioprojective",
    "regroup_any_array",
]


def only_numbers(string):
    """Remove all non-numeric characters from a string.

    For example, \"1a2H 3;\" becomes \"123\".
    """
    return re.sub(r"\D", "", string)


def regroup_any_array(
    data: np.ndarray | u.Quantity,
    old_bins: np.ndarray | u.Quantity,
    new_bins: np.ndarray | u.Quantity,
    combine_by: str | None = None,
):
    """Takes any array of data in old_bins space and rebins along data
    array axis==0 to have new_bins.

    Can specify how the bins are combined: [\"sum\", \"mean\", \"quadrature\"].

    Parameters
    ----------
    data, old_bins, new_bins : `~numpy.ndarray` | `~astropy.units.Quantity`
            Array of the data, current bins for the data, and new bins
            for the data. Shape of the bin arrays should be `(N, 2)`
            where `N` is the length of `data`.

    combine_by : string
            Defines how to combine multiple bins along axis 0. E.g., \"sum\"
            adds the data, \"mean\" averages the data, and \"quadrature\"
            sums the data in quadrature. If `None`, then \"sum\" is used.
            Default: None

    Returns
    -------
    The grouped data array.
    """
    # sanatise inputs
    combine_by = "sum" if combine_by is None else combine_by
    data, du = _get_val_and_unit(data)
    old_bins, obu = _get_val_and_unit(old_bins)
    new_bins, nbu = _get_val_and_unit(new_bins)
    old_bins = _convert_old_value_to_new_unit_values(old_bins, obu, nbu)

    new_binned_data = []
    for nb in new_bins:
        # just loop through new bins and bin data from between new_bin_lower<=old_bin_lowers and old_bin_highers<new_bin_higher
        if combine_by == "sum":
            new_binned_data.append(
                np.sum(
                    data[
                        np.nonzero(
                            (nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1])
                        )
                    ],
                    axis=0,
                )
            )
        elif combine_by == "mean":
            new_binned_data.append(
                np.mean(
                    data[
                        np.nonzero(
                            (nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1])
                        )
                    ],
                    axis=0,
                )
            )
        elif combine_by == "quadrature":
            new_binned_data.append(
                np.sqrt(
                    np.sum(
                        data[
                            np.nonzero(
                                (nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1])
                            )
                        ]
                        ** 2,
                        axis=0,
                    )
                )
            )
    return np.array(new_binned_data) if du is None else np.array(new_binned_data) << du


def _get_val_and_unit(value: np.ndarray | u.Quantity | float):
    """Return the value of an object and unit if possible.

    Returns value and `None` if no unit.
    """
    return (value.value, value.unit) if isinstance(value, u.Quantity) else (value, None)


def _convert_old_value_to_new_unit_values(
    old_value: np.ndarray | float,
    old_unit: u.core.PrefixUnit | u.core.CompositeUnit,
    new_unit: u.core.PrefixUnit | u.core.CompositeUnit,
):
    """Convert the old value to new units, return old value if no units.

    Returns unitless value.
    """
    if (old_unit is not None) and (new_unit is not None):
        return ((old_value << old_unit) << new_unit).value
    return old_value


def assign_plot_settings(nustar_map):
    """Function to set plot settings that are multiple lines."""
    nustar_map.plot_settings["norm"] = IMAGE_SCALE
    nustar_map.plot_settings["cmap"] = IMAGE_COLORMAP
    return nustar_map


def pixel_coords_to_solar_coords(pixel_coords):
    """Convert solar pixel coordinates to solar heliographic coordinates."""
    pix_map_centre = 1500
    sol_x_pix = (pixel_coords - pix_map_centre) << u.pixel
    return sol_x_pix * IMAGE_RES


def arcsec_points_to_shape(*args):
    """A function to turn a list of (x, y) points into a polygon.

    The coordinates will be the vertices of the desired shape and passed
    to ``shapely.geometry.polygon.Polygon``.

    Parameters
    ----------
    *args :
        The vertices of the polygon and should be in or convertible to
        ``~atropy.units.Quantity``.
    """
    # make sure arcsecond
    # need to index [0] as it weirdly becomes a list in a list for some reason
    unit_check_list = next(a << u.arcsec for a in args)
    unitless_for_poly = [ucl.value for ucl in unit_check_list]
    return Polygon(unitless_for_poly)


def radec_to_helioprojective(ra, dec, obstime, **kwargs):
    """Convert RA-Dec to solar helioprojective coordinates for NuSTAR.

    Parameters
    ----------
    ra : `astropy.units.Quantity`
        The right ascension position in, or convertible to, degrees.

    y : `astropy.units.Quantity`
        The declination position in, or convertible to, degrees.

    time : `astropy.time.Time`
        The time corresponding to the (ra,dec) position.

    **kwargs : Passed to ``nustar_pysolar.convert._delta_solar_skyfield``

    Returns
    -------
    `tuple[astropy.units.Quantity]` :
        The solar-x and -y in arcseconds.
    """
    from nustar_pysolar.convert import _delta_solar_skyfield

    ra <<= u.deg
    dec <<= u.deg
    sol_x, sol_y = _delta_solar_skyfield(
        np.array([ra.value]) << ra.unit,
        np.array([dec.value]) << dec.unit,
        [obstime],
        **kwargs,
    )
    return sol_x[0] << u.arcsec, sol_y[0] << u.arcsec


def helioprojective_to_radec(x, y, time, **kwargs):
    """Convert solar helioprojective to RA-Dec coordinates for NuSTAR.

    This function takes the place of ``get_skyfield_position`` and
    ``get_sky_position`` from the ``nustar_pysolar.planning`` package
    and module. No matter what I did, I could not get the correct values
    from these functions so I (Kris) wrote the exact reverse of
    ``nustar_pysolar.convert._delta_solar_skyfield`` which is used in
    ``jdhsrpy.utils.radec_to_helioprojective``.

    This was tested by confirming
    ``nustar_pysolar.convert._delta_solar_skyfield`` correctly converted
    RA-Dec into helioprojective using NuSTAR maps, SAOImageDS9, and
    region REG files. Once RA-Dec -> Solar-X/Y was confirmed, these
    coordinates were converted back into RA-Dec to compare to the
    original SAOImageDS9 REG file coordinates.

    Parameters
    ----------
    x : `astropy.units.Quantity`
        The solar-x position in, or convertible to, arcseconds.

    y : `astropy.units.Quantity`
        The solar-y position in, or convertible to, arcseconds.

    time : `astropy.time.Time`
        The time corresponding to the (x,y) position.

    **kwargs : Only used to pass ``load_path`` to
        ``nustar_pysolar.utils.skyfield_ephem``

    Returns
    -------
    `tuple[astropy.units.Quantity]` :
        The Ra and Dec in degrees.
    """
    from nustar_pysolar.utils import skyfield_ephem

    load_path = kwargs.get("load_path", None)
    observer, the_sun, ts = skyfield_ephem(
        load_path=load_path, parallax_correction=True, utc=time
    )

    x <<= u.arcsec
    y <<= u.arcsec
    arcxy = np.array([x.value, y.value]) << u.arcsec

    tcheck = ts.from_astropy(time)
    astrometric = observer.at(tcheck).observe(the_sun)
    this_ra, this_dec, _ = astrometric.radec()

    sun_pos = np.array([this_ra.to(u.deg).value, this_dec.to(u.deg).value]) * u.deg
    sun_np = sun.P(time)

    # reverse East-West conversion
    delta_offset = arcxy * [-1, 1]

    # project off the Sun
    rotation_matrix = np.array(
        [[np.cos(sun_np), np.sin(sun_np)], [-np.sin(sun_np), np.cos(sun_np)]]
    )
    delta_offset <<= u.deg
    offset = np.dot(delta_offset, np.linalg.inv(rotation_matrix))

    position = offset + sun_pos

    return position


def deconstruct_degminsec(degminsec_str):
    """Will deconstruct and deg:min:sec string into just deg."""
    degminsec = degminsec_str.split(":")
    abs_deg_part = abs(int(degminsec[0]))
    deg = (
        (abs_deg_part << u.deg)
        + (int(degminsec[1]) << u.arcmin)
        + (float(degminsec[2]) << u.arcsec)
    )
    return deg * int(degminsec[0]) / abs_deg_part


def deconstruct_hrminsec(hrminsec_str):
    """Will deconstruct and hr:min:sec string into just deg."""
    return deconstruct_degminsec(hrminsec_str) * 15
