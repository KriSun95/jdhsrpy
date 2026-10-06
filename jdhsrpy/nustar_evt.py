import warnings

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import ntpath
import numpy as np
import re
import sunpy.map

from jdhsrpy import NUSTAR_EPOCH, PIXELUNIT, IMAGE_RES
from jdhsrpy.list_filters import by_good_pix, by_energy, by_gradezero
from jdhsrpy.utils import regroup_any_array, assign_plot_settings

__all__ = ["NustarEvt", "sunpos_evt", "livetime_array", "chu_state_array", "nustar_time_from_utc", "utc_from_nustar_time", "make_sunpy_map"]


class NustarEvt():
    """A class to load in and work with NuSTAR EVT and EVT related files.

    For full functionality, this class expects an EVT file with 
    coordinates already covnerted to solar X/Y. However, this only 
    affects the image plotting methods.

    See ``~jdhsrpy.nustar_evt.sunpos_evt`` for more details on obtaining 
    the necessary file.
    
    Parameters
    ----------
    evt_filename : `str`
        The path and file name to the EVT file.
    """

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
        """Filter the event list.
        
        Filter by energy (default 2.5<=E<80), only grade 0, 
        and remove events from bad pixels from the given FPM."""
        fpm = self.fpm if fpm is None else fpm
        energy_low = 2.5 if energy_low is None else energy_low
        energy_high = 80 if energy_high is None else energy_high
        _zero = by_gradezero(event_data)
        _eng = by_energy(_zero, energy_low=energy_low, energy_high=energy_high)
        return by_good_pix(_eng, fpm=fpm)

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
        """Hacky fix for nustar_pysolar not updating the header units with the value."""
        ############*********** this is a hacky fix but will do for now ***********############
        # if Python code is used for the sunpos file creation the re-written header keywords might not save properly, so...
        if not np.allclose([self.evt_header['TCDLT13'], self.evt_header['TCDLT14']], [2.5, 2.5], atol=1e-1):
            self.evt_header['TCDLT13'] = IMAGE_RES # x
            self.evt_header['TCDLT14'] = IMAGE_RES # y

    def count_time_profile_array(self, event_data=None, time_bins=None, time_binning=None, start_time=None, end_time=None):
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
            return self.count_time_profile_array_from_start_stop_binning(event_data=event_data, time_binning=time_binning, start_time=start_time, end_time=end_time)
        return self.count_time_profile_array_from_binning_array(event_data=event_data, time_bins=time_bins)

    def count_time_profile_array_from_start_stop_binning(self, event_data=None, time_binning=None, start_time=None, end_time=None):
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
        return utc_from_nustar_time(time_bins<<u.second), counts<<u.ct
    
    def count_time_profile_array_from_binning_array(self, event_data=None, time_bins=None):
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
        counts, time_bins = np.histogram(event_data['TIME'], time_bins)
        return utc_from_nustar_time(time_bins<<u.second), counts<<u.ct

    def rate_time_profile_array(self, hk_filename, event_data=None, time_bins=None, time_binning=None, start_time=None, end_time=None):
        """Get counts and time bins for plotting a NuSTAR time profile.
        
        Parameters
        ----------
        hk_filename : `str`
            The FITS file containing the livetime information. Should be of 
            the form `nu<OBSID><FPM>_fpm.hk` and may be in the HK folder.
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
        time_bins, counts = self.count_time_profile_array(event_data=event_data, 
                                                          time_bins=time_bins, 
                                                          time_binning=time_binning, 
                                                          start_time=start_time, 
                                                          end_time=end_time)
        lvt_times, hk_livetimes = livetime_array(hk_filename)

        format_c_times = np.hstack((time_bins[:-1][:,None], time_bins[1:][:,None]))
        format_lvt_times = np.hstack((lvt_times[:-1][:,None], lvt_times[1:][:,None]))
        new_livetimes = regroup_any_array(data=hk_livetimes, 
                                          old_bins=format_lvt_times, 
                                          new_bins=format_c_times, 
                                          combine_by="mean")
        time_diff = (time_bins[1:] - time_bins[:-1]).sec << u.s
        return time_bins, counts/(new_livetimes.value*time_diff)

    def _native_map(self, event_data=None):
        """Return a Sunpy map of the NuSTAR data."""
        event_data = self.cleaned_evt_data if event_data is None else event_data
        header = self.evt_header
        start_t = str(utc_from_nustar_time(event_data["TIME"][0]<<u.second))
        end_t = str(utc_from_nustar_time(event_data["TIME"][-1]<<u.second))
        return make_sunpy_map(event_data, header, start_t=start_t, end_t=end_t)

    def full_disk_bounds(self):
        """Obtain the on-disk bounds of the NuSTAR event list."""
        bottom_left = [-1200, -1200] << u.Unit("arcsec")
        top_right = [1200, 1200] << u.Unit("arcsec")
        return bottom_left, top_right

    def full_disk_map(self, event_data=None):
        """Return a Sunpy map of the NuSTAR data."""
        nustar_map = self._native_map(event_data=event_data)
        bottom_left, top_right = self.full_disk_bounds()
        return get_submap(nustar_map, bottom_left, top_right)

    def field_of_view_bounds(self, sunpy_map_obj):
        """Obtain the on-disk bounds of the NuSTAR event list."""
        pixel_pos = np.argwhere(sunpy_map_obj.data>0) * u.pixel
        hpc_max = sunpy_map_obj.wcs.pixel_to_world(pixel_pos[:, 1], pixel_pos[:, 0])
        bottom_left = [np.min(hpc_max.Tx.value), 
                       np.min(hpc_max.Ty.value)] << u.Unit(PIXELUNIT)
        top_right = [np.max(hpc_max.Tx.value), 
                     np.max(hpc_max.Ty.value)] << u.Unit(PIXELUNIT)
        return bottom_left, top_right

    def field_of_view_map(self, event_data=None):
        """Return a map that is bound by where the event counts are."""
        nustar_map = self._native_map(event_data=event_data)
        bottom_left, top_right = self.field_of_view_bounds(nustar_map)
        return get_submap(nustar_map, bottom_left, top_right)

    def nustar_map(self, bottom_left, top_right, event_data=None):
        """Create a map with custom bounds.
        
        Parameters
        ----------
        bottom_left : `~astropy.units.Quantity`
            The bottom left coordinates of the submap. Unit aware (e.g., 
            arc-seconds)

        top_right : `~astropy.units.Quantity`
            The top right coordinates of the submap. Unit aware (e.g., 
            arc-seconds)
        """
        nustar_map = self._native_map(event_data=event_data)
        return get_submap(nustar_map, bottom_left, top_right)

def sunpos_evt(file, load_path=None):
    """Convert a .evt NuSTAR file to a _sunpos.evt file.
    
    The new file will contain Solar-X/-Y coordinates for each event.
    """
    # importing this causes the docs to fail for this file
    import nustar_pysolar
    load_path = "./" if load_path is None else load_path
    nustar_pysolar.convert.convert_file(file, load_path=load_path)

def chu_state_array(chu_filename):
    """Get an array of camera head unit states for NuSTAR from a file.

    [1] https://github.com/ianan/nustar_sac/blob/master/idl/load_nschu.pro
    [2] https://github.com/NuSTAR/nustar_solar/blob/master/depricated/solar_mosaic_20150429/read_chus.pro

    Parameters
    ----------
    chu_filename : `str`
        The FITS file containing the CHU information. Should be of the 
        form `nu<OBSID>_chu123.fits` and may be in the HK folder.

    Returns
    -------
    The times, CHU states and correesponding labels.
    """
    #not self.chu_filename as fits.open needs to know the full path to the file
    with fits.open(chu_filename) as hdulist:
        data1 = hdulist[1].data
        data2 = hdulist[2].data
        data3 = hdulist[3].data

    # easier to work with numpy arrays later
    data_c1 = np.array(data1)
    data_c2 = np.array(data2)
    data_c3 = np.array(data3)

    maxres = 20
    
    for chu_num, dat in enumerate([data_c1, data_c2, data_c3]):
        chu_bool = ((dat['VALID']==1) & 
                    (dat['RESIDUAL']<maxres) &
                    (dat['STARSFAIL']<dat['OBJECTS']) &
                    (dat['CHUQ'][:,3]!=1))
        chu_01 = chu_bool*1 # change true/false into 1/0

        chu_mask = chu_01* (chu_num+1)**2 # give each chu a unique number that when it is added to another it gives a unique chu combo, like file permissions

        if chu_num == 0:
            chu_all = chu_mask # after chu 1 file have an array with 1s and 0s
        else:
            chu_all += chu_mask # after the others (chu2 and chu3) have an array with 1,4,9,5,10,13,14
    
    # last data array in the for loop can give the time, no. of seconds from 1-Jan-2010
    chu_time = dat['TIME']

    tick_labels = ['1', '2', '12', '3', '13', '23', '123'] 
    tick_values = [100, 101, 102, 103, 104, 105, 106]
    tick_label_map = dict(zip(tick_labels, tick_values))

    # reassigned values are at 100, etc. as to not accidently double sort the values again
    # e.g. if mask value was changed to 10, then if it was accidently run again it would get sorted into chu state 13 etc.
    chu_all[chu_all == 1] = tick_label_map["1"] #chu1 # mask value in array is changed to chu state, e.g. mask value=5, chu state is 12, and value 102
    chu_all[chu_all == 4] = tick_label_map["2"] #chu2 
    chu_all[chu_all == 5] = tick_label_map["12"] #chu12
    chu_all[chu_all == 9] = tick_label_map["3"] #chu3
    chu_all[chu_all == 10] = tick_label_map["13"] #chu13
    chu_all[chu_all == 13] = tick_label_map["23"] #chu23
    chu_all[chu_all == 14] = tick_label_map["123"] #chu123

    chu_time = chu_time[chu_all > 0] # if there is still no chu assignment for that time then remove
    chu_all = chu_all[chu_all > 0]

    chu_times = utc_from_nustar_time(chu_time<<u.s) 

    return Time(chu_times), chu_all, tick_label_map

def livetime_array(hk_filename):
    """Get an array of the livetime fractions for NuSTAR from a file.

    Parameters
    ----------
    hk_filename : `str`
        The FITS file containing the livetime information. Should be of 
        the form `nu<OBSID><FPM>_fpm.hk` and may be in the HK folder.

    Returns
    -------
    The times and livetimes.
    """
    with fits.open(hk_filename) as hdulist:
        hk_data = hdulist[1].data

    lvt_times = utc_from_nustar_time(hk_data['time']<<u.s) 
    hk_livetimes = hk_data['livetime']<<u.percent

    return Time(lvt_times), hk_livetimes

def nustar_time_from_utc(utc_time:Time):
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
    return (utc_time - NUSTAR_EPOCH).sec

def utc_from_nustar_time(nustar_time:u.Quantity):
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
    return NUSTAR_EPOCH + (nustar_time << u.s)

def make_sunpy_map(evtdata, hdr, **kwargs):
    """ Make a sunpy map based on the NuSTAR data.
    
    Parameters
    ----------
    evtdata: FITS data structure
        This should be an hdu.data structure from a NuSTAR FITS file.

    hdr: FITS header containing the astrometric information

    **kwargs :
        Added to the header.
    """

    # Parse Header keywords
    for field in hdr.keys():
        if field.find('TYPE') != -1:
            if hdr[field] == 'X':
                xval = field[5:8]
            if hdr[field] == 'Y':
                yval = field[5:8]
        
    min_x = hdr['TLMIN'+xval]
    min_y = hdr['TLMIN'+yval]
    max_x = hdr['TLMAX'+xval]
    max_y = hdr['TLMAX'+yval]

    delx = abs(hdr['TCDLT'+xval])

    x = evtdata['X'][:]
    y = evtdata['Y'][:]
    met = evtdata['TIME'][:]*u.s
    mjdref = hdr['MJDREFI']

    mid_obs_time = Time(mjdref*u.d+met.mean(), format = 'mjd')
    sta_obs_time = Time(mjdref*u.d+met.min(), format = 'mjd')

    # Add in the exposure time (or livetime), just a number not units of seconds 
    exp_time = hdr['EXPOSURE']

    # Assume X and Y are the same size
    resample = 1.0
    scale = delx * resample
    bins = int((max_x - min_x) / (resample))

    H, _, _ = np.histogram2d(y, x, bins=bins, range = [[min_y,max_y], [min_x, max_x]])

    #Normalise the data with the exposure (or live) time?
    pixluname='count'

    dict_header = {
        "DATE-OBS": sta_obs_time.iso, 
        "EXPTIME": exp_time,
        "CDELT1": scale,
        "NAXIS1": bins,
        "CRVAL1": 0.,
        "CRPIX1": bins*0.5,
        "CUNIT1": PIXELUNIT,
        "CTYPE1": "HPLN-TAN",
        "CDELT2": scale,
        "NAXIS2": bins,
        "CRVAL2": 0.,
        "CRPIX2": bins*0.5 + 0.5,
        "CUNIT2": PIXELUNIT,
        "CTYPE2": "HPLT-TAN",
        "PIXLUNIT": pixluname,
        "HGLT_OBS": sunpy.coordinates.sun.B0(mid_obs_time), #get_sun_B0(mid_obs_time),
        "HGLN_OBS": 0,
        "RSUN_OBS": sunpy.coordinates.sun.angular_radius(mid_obs_time).value, #sun.solar_semidiameter_angular_size(mid_obs_time).value,
        "RSUN_REF": sunpy.sun.constants.radius.value, #sun.constants.radius.value,
        # Assumes dsun_obs in m if don't specify the units, so give units
        "DSUN_OBS": sunpy.coordinates.sun.earth_distance(mid_obs_time).value*u.astrophys.au, #get_sunearth_distance(mid_obs_time).value*u.astrophys.au
    }
    dict_header |= kwargs
    # For some reason the DSUN_OBS crashed the save...
    
    header = sunpy.util.MetaDict(dict_header)

    nustar_map = sunpy.map.Map(H, header)
    nustar_map = assign_plot_settings(nustar_map)
    
    return nustar_map

def get_submap(sunpy_map_obj, bottom_left:u.Quantity, top_right:u.Quantity):
    """Create and return a submap of a Sunpy map object.
    
    Parameters
    ----------
    sunpy_map_obj : ``
        The map object.
    
    bottom_left : `~astropy.units.Quantity`
        The bottom left coordinates of the submap. Unit aware (e.g., 
        arc-seconds)
    
    top_right : `~astropy.units.Quantity`
        The top right coordinates of the submap. Unit aware (e.g., 
        arc-seconds)
    """
    # fix really large plot, instead of going from -3600 to 3600 in x and y
    bl = SkyCoord(*bottom_left, frame=sunpy_map_obj.coordinate_frame)
    tr = SkyCoord(*top_right, frame=sunpy_map_obj.coordinate_frame)
    _submap = sunpy_map_obj.submap(bl,top_right=tr) # submaps act differently to normal maps
    return assign_plot_settings(sunpy.map.Map(_submap.data, _submap.meta))

def normalize_map(sunpy_map_obj, hk_filename):
    """Normalize the given map data with its EXPOSURE header keyword.
    
    Updates the PIXLUNIT keyword value too.
    """
    if sunpy_map_obj.meta["PIXLUNIT"].endswith("/s"):
        warnings.warn("PIXLUNIT in map object already ends in `/s`. Returning original map.")
        return sunpy_map_obj
    
    lvt_exp = get_livetime_exposure(Time(sunpy_map_obj.meta["start_t"], format='isot',scale='utc'), 
                                    Time(sunpy_map_obj.meta["end_t"], format='isot',scale='utc'), 
                                    hk_filename)

    new_data = sunpy_map_obj.data/lvt_exp.value
    sunpy_map_obj.meta["PIXLUNIT"] = f"{sunpy_map_obj.meta["PIXLUNIT"]}/s"
    return assign_plot_settings(sunpy.map.Map(new_data, sunpy_map_obj.meta))

def get_livetime_exposure(start, end, hk_filename):
    """Use the HK file and time span to get the livetime corrected exposure."""
    lvt_times, hk_livetimes = livetime_array(hk_filename)
    new_time_bins = [[start, end]]

    format_lvt_times = np.hstack((lvt_times[:-1][:,None], lvt_times[1:][:,None]))
    new_livetimes = regroup_any_array(data=hk_livetimes, 
                                        old_bins=format_lvt_times, 
                                        new_bins=new_time_bins, 
                                        combine_by="mean")
    time_diff = (end - start).sec << u.s
    return (new_livetimes.value*time_diff)

def draw_grid(sunpy_map_obj, axes):
    """Draw a grid representing the Sun."""
    sunpy_map_obj.draw_limb(axes=axes, color='black',linewidth=1,linestyle='dashed', zorder=0)
    sunpy_map_obj.draw_grid(axes=axes, color='grey', alpha=0.5,linewidth=1,linestyle='dashed', zorder=0, annotate=False)
