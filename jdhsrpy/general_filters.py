__all__ = [
    "by_value_range",
]


def by_value_range(array, start, stop, array_to_filter_against=None):
    """Do some basic filtering on an array with a start and stop.

    Parameters
    ----------
    array: array-like
        This should be an hdu.data structure from a NuSTAR FITS file.

    start : value-like
        The start of the range to slice ``array`` (inclusive).

    stop : value-like
        The stop of the range to slice ``array`` (exclusive).

    array_to_filter_against : array-like
        If not None, then indices from this array within ``start`` and
        ``stop`` will be used to select values in ``array``.

    Returns
    -------
    The sliced list.
    """
    farray = array if array_to_filter_against is None else array_to_filter_against
    gen_filter = (farray >= start) & (farray < stop)
    inds = (gen_filter).nonzero()
    goodinds = inds[0]

    return array[goodinds]
