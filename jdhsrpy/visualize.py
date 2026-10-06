from astropy.time import Time
from astropy.visualization import time_support
from datetime import datetime
import matplotlib.pyplot as plt

__all__ = ["time_profile_plot", "vertical_line_of_time"]

def time_profile_plot(times, ydata, axes=None, **kwargs):
    """Plot a time profile (on a given axes if given).
    
    Parameters
    ----------
    times : `~astropy.time.Time`|`~datetime.datetime`
        The time edges.

    ydata : `list`|`~numpy.ndarray`
        The data varying with time. Length will be one less the `times`. 
        length.

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    time_support(format='unix_tai')
    axes.stairs(ydata, to_plotting_datetimes(times), **kwargs)
    return axes

def to_plotting_datetimes(times):
    """Convert times to ones that plot well."""
    if isinstance(times, Time):
        return Time(times, format='unix_tai',scale='utc').datetime
    elif isinstance(times, datetime):
        return times
    else:
        raise ValueError("The `times` variable is not supported.")

def vertical_line_of_time(time:str|Time|datetime, axes=None, **kwargs):
    """Use `axvline` to plot `time` input.
    
    Parameters
    ----------
    time : `str`|`astropy.time.core.Time`|`datetime`
        If `str` then need ISOT format. E.g., 2021-11-20T02:25:30.

    Returns
    -------
    Axes input.
    """
    axes = plt if axes is None else axes
    if isinstance(time, str):
        converted_time = Time(time, format='isot',scale='utc').datetime
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

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    time_support(format='unix_tai')
    axes.plot(to_plotting_datetimes(times), ydata, "x", **kwargs)
    axes.axes.set_yticks(list(label_map.values()))
    axes.axes.set_yticklabels(list(label_map.keys()))
    return axes

def spectrum_plot(x, y, xerr=None, yerr=None, axes=None, **kwargs):
    """Functoin to plot a spectrum."""
    axes = plt.gca() if axes is None else axes
    _defaults = {"color":"k", 
                 "fmt":".", 
                 "markersize":0.01}
    _plot_settings = _defaults | kwargs
    axes.errorbar(x, y, xerr=xerr, yerr=yerr, **_plot_settings)
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

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    _defaults = {"drawstyle":"steps-mid",
                 }
    _plot_settings = _defaults | kwargs
    plt.semilogy(to_plotting_datetimes(times), ydata, **_plot_settings)
    return axes

def spectrum_plot(energy_bins, ydata, ydata_err, axes=None, **kwargs):
    """Plot a time profile (on a given axes if given).
    
    Parameters
    ----------
    energy_bins : `~astropy.unit.Quantity`
        The energy bin edges. Should be one more than the `ydata` legnth.

    ydata :  `~astropy.unit.Quantity`
        The spectrum.

    Returns
    -------
    Axes input.
    """
    axes = plt.gca() if axes is None else axes
    _defaults = {"color":"k",
                 "fmt":".",
                 "markersize":0.01,
                 }
    _plot_settings = _defaults | kwargs

    energy_mids = (energy_bins[:-1] + energy_bins[1:])/2
    energy_half_width = (energy_bins[1:] - energy_bins[:-1])/2

    ## plot spectrum created from Python
    axes.errorbar(energy_mids, ydata, xerr=energy_half_width, yerr=ydata_err, **_plot_settings)

    axes.set_ylabel(f"{ydata.unit:latex}")
    axes.tick_params(axis='both', which='both', direction='in')
    axes.set_yscale('log')
    axes.set_xlabel(f"Energy [{energy_bins.unit:latex}]")
    return axes
