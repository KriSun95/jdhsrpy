import re

import astropy.units as u
import numpy as np

from jdhsrpy import IMAGE_COLORMAP, IMAGE_SCALE

__all__ = ["only_numbers", "regroup_any_array", "assign_plot_settings"]

def only_numbers(string):
    """Remove all non-numeric characters from a string.
    
    For example, \"1a2H 3;\" becomes \"123\". 
    """
    return re.sub(r"\D", "", string)

def regroup_any_array(data:np.ndarray|u.Quantity, old_bins:np.ndarray|u.Quantity, new_bins:np.ndarray|u.Quantity, combine_by:str|None=None):
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
    new_bins, nbu  = _get_val_and_unit(new_bins)
    old_bins = _convert_old_value_to_new_unit_values(old_bins, obu, nbu)

    new_binned_data = []
    for nb in new_bins:
        # just loop through new bins and bin data from between new_bin_lower<=old_bin_lowers and old_bin_highers<new_bin_higher
        if combine_by == "sum":
            new_binned_data.append(
                np.sum(data[np.nonzero((nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1]))], axis=0)
            )
        elif combine_by == "mean":
            new_binned_data.append(
                np.mean(data[np.nonzero((nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1]))], axis=0)
            )
        elif combine_by == "quadrature":
            new_binned_data.append(
                np.sqrt(np.sum(data[np.nonzero((nb[0] <= old_bins[:, 0]) & (nb[-1] >= old_bins[:, -1]))] ** 2, axis=0))
            )
    return np.array(new_binned_data) if du is None else np.array(new_binned_data) << du

def _get_val_and_unit(value:np.ndarray|u.Quantity|float|int):
    """Return the value of an object and unit if possible.

    Returns value and `None` if no unit.
    """
    return (value.value, value.unit) if isinstance(value, u.Quantity) else (value, None)

def _convert_old_value_to_new_unit_values(old_value:np.ndarray|float|int, old_unit:u.core.PrefixUnit|u.core.CompositeUnit, new_unit:u.core.PrefixUnit|u.core.CompositeUnit):
    """Convert the old value to new units, return old value if no units.

    Returns unitless value.
    """
    if (old_unit is not None) and (new_unit is not None):
        return ((old_value<<old_unit)<<new_unit).value
    return old_value

def assign_plot_settings(nustar_map):
    """Function to set plot settings that are multiple lines."""
    nustar_map.plot_settings['norm'] = IMAGE_SCALE
    nustar_map.plot_settings['cmap'] = IMAGE_COLORMAP
    return nustar_map