"""
File to host code that downloads SDO data.
"""

import logging
import os

import astropy.units as u
from sunpy.net import Fido
from sunpy.net import attrs as a


def sdo_download(start_time, end_time, directory=None, wave_values=None, get_hmi=False):
    """Allows the download of SDO/AIA and /HMI data.

    If the files already exist in the same directory, they will not be
    downloaded again.

    *** Due to downloading stuff, this function might need to be run
    several times. ***

    Parameters
    ----------
    start_time : `str`, `None`
        The start time (UTC) for the data. In the following format; e.g.,
        "2024-04-17T21:57:00".

    end_time : `str`, `None`
        The end time (UTC) for the data. In the following format; e.g.,
        "2024-04-17T22:30:00".

    directory : `str`
        The directory where the SDO files will be saved. Default is the
        current directory. Default is `os.getcwd()` if `None`.
        Default: None

    wave_values : `list[int]`, `None`
        A list of the wavelengths (in angstroms) to be downloaded. For
        example, [94, 131, 171, 211, 1600, 1700] which is the value
        taken if input is `None`.
        Default: None

    get_hmi : `Bool`
        A boolean flag where, if True, HMI LOS data will also be
        downloaded.
        Default: False

    Returns
    -------
    : `list[str]`
        A list object containing the names of the files that failed to
        download.

    Examples
    --------
    # Default behaviour for the FOXSI M1 flare time
    >>> sdo_download()

    # Extend the time interval, download a subset of wavelengths,
    # define a specific directory, and get the SDO/HMI data as well
    >>> sdo_download(start_time="2024-04-17T21:30:00", end_time="2024-04-17T22:45:00", directory="path/to/save/directory/", wave_values=[94,171], get_hmi=True)

    # Only download SDO/HMI LOS data
    >>> sdo_download(..., wave_values=[], get_hmi=True)
    """
    # set-up defaults
    directory = os.getcwd() if directory is None else directory
    wave_values = (
        [94, 131, 171, 193, 211, 304, 335, 1600, 1700]
        if wave_values is None
        else wave_values
    )

    # get ready for Sunpy
    time = a.Time(start_time, end_time)
    waves = [a.Wavelength(wv) for wv in wave_values << u.angstrom]
    wave_dirs = [f"{wv}angstrom" for wv in wave_values]
    instruments = [a.Instrument("AIA")] * len(waves)

    # check if HMI is wanted
    if get_hmi:
        waves.append(a.Physobs("LOS_magnetic_field"))
        wave_dirs.append("hmi")
        instruments.append(a.Instrument("hmi"))

    needed_files = []
    # for aia (and hmi?)
    for inst, wave, wd in zip(instruments, waves, wave_dirs):
        curr_wdir = os.path.join(directory, wd)
        logging.info(f"Doing {curr_wdir}")
        os.makedirs(curr_wdir, exist_ok=True)  # make the directory if it isn't there

        search = (time, inst, wave)
        fetch = {"path": os.path.join(curr_wdir, "{file}")}

        filepaths = fido_download(search, fetch)

        # well try again with failed files
        needed_files = handle_retries(filepaths, needed_files=needed_files)

    logging.info("Files Needed are:")
    logging.info(needed_files)

    return needed_files


def fido_download(search, fetch=None):
    """Handle the common download elements for `sunpy.Fido`.

    Parameters
    ----------
    search : `tuple`
        The inputs given to the `Fido.search` query.

    fetch : `dict`
        Kyeword arguments to be passed to the `Fido.fetch` method. If
        `None` then an empty `dict` object is passed.
        Default: None

    Returns
    -------
    : `sunpy.UnifiedResponse` object
        A list object containing the names of the files that failed to
        download.

    Examples
    --------
    ## Download some SDO/AIA data
    >>> from sunpy.net import attrs as a
    >>> search = (a.Time("2024-04-17T21:30:00",
                         "2024-04-17T22:30:00"),
                  a.Instrument('AIA'),
                  a.Wavelength(wave))
    >>> fetch = {"path":"path/to/save/{file}"}
    >>> filepaths = fido_download(search, fetch)

    ## Download some SDO/HMI data
    >>> from sunpy.net import attrs as a
    >>> search = (a.Time("2024-04-17T21:30:00",
                         "2024-04-17T22:30:00"),
                  a.Instrument('hmi'),
                  a.Physobs('LOS_magnetic_field'))
    >>> fetch = {"path":"path/to/save/{file}"}
    >>> filepaths = fido_download(search, fetch)
    """
    unifresp = Fido.search(*search)

    fetch = {} if fetch is None else fetch
    filepaths = Fido.fetch(unifresp, **fetch)

    return filepaths


def handle_retries(filepaths, tries=5, needed_files=None):
    """If a file(s) isn't downloaded straight away, retry the download.

    Parameters
    ----------
    filepaths : `sunpy.UnifiedResponse` object
        An object used to store the failed downloaded files.

    tries : `int`
        The number of times to retry downloading the given files.
        Default: 5

    needed_files : `list[str]`
        A list of files already needed to be returned if failed after
        retrying `tries` times.For example, [] which is the value taken
        if input is `None`. If this is a list object then the failed
        files will be appended to this list.
        Default: None

    Returns
    -------
    : `list[str]`
        A list object containing the file names of the failed downloads
        after retrying.

    Examples
    --------
    ## Retry failed downloads and return a list still un-downloaded
    >>> needed_files = handle_retries(filepaths)

    ## If a list of un-downloaded already exists then pass to append
    >>> needed_files = handle_retries(filepaths, needed_files=needed_files)
    """

    needed_files = [] if needed_files is None else needed_files

    for _ in range(tries):
        logging.info(f"Trying again for {filepaths.errors}")
        filepaths = Fido.fetch(filepaths)
        if filepaths.errors == []:
            break

    if filepaths.errors != []:
        needed_files.append(filepaths.errors)

    return needed_files
