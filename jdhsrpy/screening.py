from astropy.io import fits
from astropy.time import Time
import ntpath
import os

from jdhsrpy import NUSTAR_EPOCH
from jdhsrpy.filter import in_time_range_inds

__all__ = ["time_filtered_evt_file", "make_gti_file"]

def time_filtered_evt_file(evt_file, time_range=None, save_dir=None, **kwargs):
    """Takes a .evt file and filters the events list to a given time range. 
    Only for region selection, do not use directly with spectral fitting 
    software.
    
    Parameters
    ----------
    file : Str
            File (or directory/file) of the .evt file to be filtered by time.
    
    time_range : list
            A list of length 2 with the start and end date and time. Must 
            be given in a specific format, e.g. time_range=['2018/09/10, 16:22:30', '2018/09/10, 16:24:30'].
            Default: None
            
    save_dir : Str
            String of the directory for the filtered file to be saved.
            Default: None
            
    Returns
    -------
    Creates a new file file with '_time_filtered' before the file extension 
    and returns the name of the new file.
    """
    
    if time_range is None:
        print('No time_range given. Nothing will be done.')
        return
    
    directory, filename = ntpath.split(evt_file)
    tf_name = f"{filename[:-4]}_time_filtered.evt"
    save_dir = directory if save_dir is None else save_dir
    new_file_name = os.path.join(save_dir, tf_name) 
    
    with fits.open(evt_file) as hdulist:
        evtdata = hdulist[1].data # data to be filtered

        timeinds = in_time_range_inds(evtdata, time_range) # picks events inside time range
        evt_in_time = evtdata[timeinds]

        hdulist[1].data = evt_in_time # replaces this hdu with the filtered events list
        hdulist.writeto(new_file_name, **kwargs) # saves the edited file, original stays as is

    return new_file_name

def make_gti_file(gti_file, save_name, good_time_interval, **kwargs):
    """Create a new GTI file with your own start and stop time.
    
    Parameters
    ----------
    gti_file : str
            The original GTI file for the NuSTAR observation. E.g., gti_file="./nu80414201001B06_gti.fits".
    
    save_name : str
            The name of your new GTI file. E.g., save_name = "./new_gti.fits"
            
    good_time_interval : [str, str]
            A list made of 2 string with the start and stop time for your good time interval. 
            E.g., good_time_interval = ['2018-09-09 09:13:36', '2018-09-09 09:17:00'] will produce 
            a new GTI file with data START='2018-09-09 09:13:36' and STOP='2018-09-09 09:17:00'.
            
    Returns
    -------
    None.

    Example
    -------
    gtiFile = "./nu80414201001B06_gti.fits"
    saveName = "./new_gti.fits"
    goodTimeInterval = ['2018-09-09T09:13:36', '2018-09-09T09:17:00']

    make_gti(gtiFile, saveName, goodTimeInterval, overwrite=True)
    """

    # custom GTI
    time1 = Time(good_time_interval[0], format='isot',scale='utc')
    time2 = Time(good_time_interval[1], format='isot',scale='utc')

    # now open the file, change the start and stop times and save to new file without changing original
    with fits.open(gti_file) as hdulist:

        # it turns out changing TSTART and TSTOP in hdulist[1].header does nothing, need to change the actual data
        hdulist[1].data["START"] = (time1 - NUSTAR_EPOCH).sec # replaces this hdu with the new time
        hdulist[1].data["STOP"] = (time2 - NUSTAR_EPOCH).sec 

        hdulist.writeto(save_name, **kwargs) # saves the edited file, original stays as is. 
        #overwrite=True overwrites save_name file if it's there, not original GTI file
