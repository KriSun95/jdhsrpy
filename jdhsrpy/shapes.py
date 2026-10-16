import astropy.units as u
from shapely.geometry import Point

from jdhsrpy.utils import arcsec_points_to_shape

__all__ = [
    "circle_shape",
    "rectangle_shape",
]


def rectangle_shape(bottom_left, top_right):
    """Return a ``shapely.geometry.polygon.Polygon`` rectangle shape.

    Parameters
    ----------
    bottom_left : `~astropy.units.Quantity`
        The bottom left corner (X, Y) of the desired box. Should be
        convertible to arcseconds.

    top_right : `~astropy.units.Quantity`
        The top right corner (X, Y) of the desired box. Should be
        convertible to arcseconds.
    """
    vertices = [
        bottom_left << u.arcsec,
        [bottom_left[0], top_right[1]] << u.arcsec,
        top_right << u.arcsec,
        [top_right[0], bottom_left[1]] << u.arcsec,
    ]
    return arcsec_points_to_shape(vertices)


def circle_shape(center, radius):
    """Return a ``shapely.geometry.polygon.Polygon`` circle shape.

    Parameters
    ----------
    center : `~astropy.units.Quantity`
        Thecenter (X, Y) of the desired circle. Should be convertible to
        arcseconds.

    radius : `~astropy.units.Quantity`
        The radius R of the desired circle. Should be convertible to
        arcseconds.
    """
    _center = center << u.arcsec
    _radius = radius << u.arcsec
    return Point(*_center.value).buffer(_radius.value)
