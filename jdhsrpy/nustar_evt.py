import warnings

from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import ntpath
import numpy as np
import re

from jdhsrpy import NUSTAR_EPOCH
from jdhsrpy.filter import bad_pix, by_energy, gradezero

__all__ = ["NUSTAR_EPOCH", "NustarEvt", "sunpos_evt"]


class NustarEvt():
    """A class to load in and work with NuSTAR EVT and EVT related files.
    
    Parameters
    ----------
    evt_filename : `str`
        The path and file name to the EVT file.
    """

    # nustar times are measured in seconds from this date
    nustar_epoch = NUSTAR_EPOCH 

    def __init__(self, evt_filename:str):
        self._check_sunpos(evt_filename)

        #extract the data within the provided parameters
        with fits.open(evt_filename) as hdulist: 
            self.evt_data = hdulist[1].data
            self.evt_header = hdulist[1].header
        self._hacky_pixel_scale_fix()

        self.directory, self.filename = ntpath.split(evt_filename)
        #search for 2 digits, a non-digit, then 2 digits again
        self.fpm = re.compile(r'\d{2}\D\d{2}').findall(self.filename)[0][2]

        self.cleaned_evt_data = self.clean_evt(self.evt_data)

    def clean_evt(self, event_data, fpm=None, energy_low=None, energy_high=None):
        fpm = self.fpm if fpm is None else fpm
        energy_low = 2.5 if energy_low is None else energy_low
        energy_high = 80 if energy_high is None else energy_high
        goodzero = gradezero(event_data)
        _zero = event_data[goodzero]
        goodeng = by_energy(_zero, energy_low=energy_low, energy_high=energy_high)
        _eng = _zero[goodeng]
        goodpix = bad_pix(_eng, fpm=fpm)
        return _eng[goodpix]

    def _check_sunpos(self, evt_filename):
        """Check for \"sunpos\" in the EVT file name."""
        # for a sunpy map object to be made then the file has to be positioned on the Sun
        if "sunpos" not in evt_filename:
            warnings.warn(
                """
                The provided file does not look like a sunpos file. Limited 
                functionality with solar coordinate related functions. See 
                ``~jdhsrpy.nustar_evt.sunpos_evt`` for help.
                """
                )

    def _hacky_pixel_scale_fix(self):
        ############*********** this is a hacky fix but will do for now ***********############
        # if Python code is used for the sunpos file creation the re-written header keywords might not save properly, so...
        if not np.allclose([self.evt_header['TCDLT13'], self.evt_header['TCDLT14']], [2.5, 2.5], atol=1e-1):
            self.evt_header['TCDLT13'] = 2.45810736 # x
            self.evt_header['TCDLT14'] = 2.45810736 # y

    def nustar_time_from_utc(self, utc_time:Time):
        """Get the number of seconds from 2010-01-01 for a given UTC time.
        
        Parameters
        ----------
        utc_time : `~astropy.time.Time`
            The UTC time Astropy object.

        Returns
        -------
        : `~astropy.units.Quantity`
            The number of seconds from the NuSTAR epoch.
        """
        return (utc_time - self.nustar_epoch).sec

    def utc_from_nustar_time(self, nustar_time:u.Quantity):
        """Get the UTC time from a number of seconds after 2010-01-01.
        
        Parameters
        ----------
        nustar_time : `~astropy.units.Quantity`
            The number of seconds from the NuSTAR epoch.

        Returns
        -------
        : `~astropy.time.Time`
            The Astropy time object in UTC.
        """
        return self.nustar_epoch + nustar_time

    def time_profile_array(self, event_data=None, time_bins=None, time_binning=None, start_time=None, end_time=None):
        """Get counts and time bins for plotting a NuSTAR time profile.
        
        Parameters
        ----------
        event_data :

            Default: None

        time_bins :
            Takes priority over `time_binning`, `start_time`, and 
            `stop_time`.
            Default: None

        time_binning :
            If `None` then a default value is used.
            Default: None

        start_time :
            If `None` then a default value is used.
            Default: None
        
        end_time :
            If `None` then a default value is used.
            Default: None
        
        Returns
        -------
        """
        if time_bins is None:
            return self.time_profile_array_from_start_stop_binning(event_data=event_data, time_binning=time_binning, start_time=start_time, end_time=end_time)
        return self.time_profile_array_from_binning_array(event_data=event_data, time_bins=time_bins)

    def time_profile_array_from_start_stop_binning(self, event_data=None, time_binning=None, start_time=None, end_time=None):
        """Get counts and time bins for plotting a NuSTAR time profile.

        User given start, stop, and/or time bin size.
        
        Parameters
        ----------
        event_data :

            Default: None

        time_binning :

            Default: None

        start_time :

            Default: None
        
        end_time :

            Default: None
        
        Returns
        -------
        """
        time_binning = 10<<u.second if time_binning is None else time_binning
        event_data = self.cleaned_evt_data if event_data is None else event_data
        start_time = np.min(event_data['TIME'])<<u.second if start_time is None else start_time
        end_time = np.max(event_data['TIME'])<<u.second if end_time is None else end_time
        start_time <<= u.second 
        end_time <<= u.second
        time_binning <<= u.second
        time_bins = np.arange(start_time.value, end_time.value+time_binning.value, time_binning.value)
        counts, time_bins =  np.histogram(event_data['TIME'], time_bins)
        return counts, self.utc_from_nustar_time(time_bins<<u.second)
    
    def time_profile_array_from_binning_array(self, event_data=None, time_bins=None):
        """Get counts and time bins for plotting a NuSTAR time profile.

        User given time bins.
        
        Parameters
        ----------
        event_data :

            Default: None

        time_bins :

            Default: None
        
        Returns
        -------
        """
        counts, time_bins =  np.histogram(event_data['TIME'], time_bins)
        return counts, self.utc_from_nustar_time(time_bins<<u.second)

def sunpos_evt(file, load_path=None):
    """Convert a .evt NuSTAR file to a _sunpos.evt file.
    
    The new file will contain Solar-X/-Y coordinates for each event.
    """
    # importing this causes the docs to fail for this file
    import nustar_pysolar
    load_path = "./" if load_path is None else load_path
    nustar_pysolar.convert.convert_file(file, load_path=load_path)
