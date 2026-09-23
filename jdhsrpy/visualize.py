from astropy.time import Time
from astropy.visualization import time_support
from datetime import datetime
import matplotlib.pyplot as plt

__all__ = ["time_profile_plot", "vertical_line_of_time"]

def time_profile_plot(times, ydata, axes=None):
    axes = plt if axes is None else axes
    if isinstance(times, Time):
        converted_times = Time(times, format='unix_tai',scale='utc').datetime
    elif isinstance(times, datetime):
        converted_times = times
    else:
        raise ValueError("The `times` variable in `time_profile_plot` is not supported.")
    time_support(format='unix_tai')
    axes.stairs(ydata, converted_times)
    return axes

def vertical_line_of_time(time:str|Time|datetime, axes=None, **kwargs):
    """Use `axvline` to plot `time` input.
    
    Parameters
    ----------
    time : `str`|`astropy.time.core.Time`|`datetime`
        If `str` then need ISOT format. E.g., 2021-11-20T02:25:30.
    """
    axes = plt if axes is None else axes
    if isinstance(time, str):
        converted_time = Time(time, format='isot',scale='utc').datetime
    elif isinstance(time, Time):
        converted_time = time.datetime
    elif isinstance(time, datetime):
        converted_time = time
    else:
        raise ValueError("The `time` variable in `vertical_line_of_time` is not supported.")
    axes.axvline(x=converted_time, **kwargs)

    