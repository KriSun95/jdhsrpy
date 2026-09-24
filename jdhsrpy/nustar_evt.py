import warnings

from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import ntpath
import numpy as np
import re

__all__ = ["NUSTAR_EPOCH", "NustarSunposEvt", "sunpos_evt", "bad_pix", "by_energy", "gradezero", "in_time_range_inds", "event_filter"]

NUSTAR_EPOCH = Time("2010-01-01T00:00:00.000", format='isot',scale='utc') 

class NustarSunposEvt():

    # nustar times are measured in seconds from this date
    nustar_epoch = NUSTAR_EPOCH 

    def __init__(self, evt_filename):
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

    def nustar_time_from_utc(self, utc_time):
        return (utc_time - self.nustar_epoch).sec

    def utc_from_nustar_time(self, nustar_time):
        return self.nustar_epoch + nustar_time

    def time_profile_array(self, event_data=None, time_binning=None, start_time=None, end_time=None):
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

def sunpos_evt(file, load_path=None):
    """Convert a .evt NuSTAR file to a _sunpos.evt file.
    
    The new file will contain Solar-X/-Y coordinates for each event.
    """
    # importing this causes the docs to fail for this file
    import nustar_pysolar
    load_path = "./" if load_path is None else load_path
    nustar_pysolar.convert.convert_file(file, load_path=load_path)

def bad_pix(evtdata, fpm):
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
    
    if fpm.find('B') == -1 :
        pix_filter = np.invert( ( (evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 16) & (evtdata['RAWY'] == 5) |
                                (evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 24) & (evtdata['RAWY'] == 22) |
                                (evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 27) & (evtdata['RAWY'] == 6) |
                                (evtdata['DET_ID'] == 2) & (evtdata['RAWX'] == 27) & (evtdata['RAWY'] == 21) |
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 22) & (evtdata['RAWY'] == 1) |
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 15) & (evtdata['RAWY'] == 3) |
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 5) & (evtdata['RAWY'] == 5) | 
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 22) & (evtdata['RAWY'] == 7) | 
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 16) & (evtdata['RAWY'] == 11) | 
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 18) & (evtdata['RAWY'] == 3) | 
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 24) & (evtdata['RAWY'] == 4) | 
                                (evtdata['DET_ID'] == 3) & (evtdata['RAWX'] == 25) & (evtdata['RAWY'] == 5) ) )
    else:
        pix_filter = np.invert( ( (evtdata['DET_ID'] == 0) & (evtdata['RAWX'] == 24) & (evtdata['RAWY'] == 24)) )


    inds = (pix_filter).nonzero()
    goodinds=inds[0]
    
    return goodinds
    
def by_energy(evtdata, energy_low=2.5, energy_high=10.):
    """ Apply energy filtering to the data.
    
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
    pi_filter = ( ( evtdata['PI']>pilow ) &  ( evtdata['PI']<pihigh))
    inds = (pi_filter).nonzero()
    goodinds=inds[0]
    
    return goodinds
    
def gradezero(evtdata):
    """ Only accept counts with GRADE==0.
        
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
    
    grade_filter = ( evtdata['GRADE'] == 0)
    inds = (grade_filter).nonzero()
    goodinds = inds[0]
    
    return goodinds

def in_time_range_inds(evtdata, tmrng):    
    """ Only include counts within a given time range.

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
    
    tstart = Time(tmrng[0], format='isot',scale='utc') 
    tend = Time(tmrng[1], format='isot',scale='utc')  
    tstart_s = (tstart - NUSTAR_EPOCH).sec #both dates are converted to number of seconds from 2010-Jan-1  
    tend_s = (tend - NUSTAR_EPOCH).sec
    tmrng = [tstart_s, tend_s] 
        
    time_filter = ( (evtdata['TIME']>tmrng[0]) & (evtdata['TIME']<tmrng[1]) )
    inds = (time_filter).nonzero()  
    goodinds=inds[0]       
 
    return goodinds 

def event_filter(evtdata, fpm='FPMA',
    energy_low=2.5, energy_high=10, tmrng = None):
    # was event_filter(evtdata, fpm='FPMA', energy_low=2.5, energy_high=10) # Kris #
    """ All in one filter module. By default applies an energy cut, 
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
    
    goodinds = in_time_range_inds(evtdata, tmrng)          
    evt_timefilter = evtdata[goodinds]     
    goodinds = bad_pix(evt_timefilter, fpm=fpm) 
    evt_badfilter = evt_timefilter[goodinds] 
    goodinds = by_energy(evt_badfilter,
                        energy_low=energy_low, energy_high = energy_high)
    evt_energy = evt_badfilter[goodinds]
    goodinds = gradezero(evt_energy)
    cleanevt = evt_energy[goodinds]
    return cleanevt
