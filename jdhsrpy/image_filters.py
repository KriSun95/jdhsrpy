
from astropy.io import fits
import numpy as np
from scipy import ndimage
from skimage import restoration
import sunpy.map

from jdhsrpy import FWHM_ARCSEC, IMAGE_RES
from jdhsrpy.utils import assign_plot_settings

__all__ = ["gaussian_filter", "deconvolve_with_array", "deconvolve_with_file"]

def gaussian_filter(sunpy_map_obj, sigma=None, **kwargs):
    """Takes a Sunpy map and applies a Guassian filter.
    
    See ``~scipy.ndimage.gaussian_filter`` for the method.

    [1] https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.gaussian_filter.html
    
    Paramteres
    ----------
    sunpy_map_obj : ``~sunpy.map.Map``
        The map for the filter to be applied to.

    sigma : scalar or sequence of scalars or None
        Standard deviation for Gaussian kernel. The standard deviations 
        of the Gaussian filter are given for each axis as a sequence, or 
        as a single number, in which case it is equal for all axes. If 
        None, then NuSTAR's FWHM of 18" and the default 2.45810736 arc/pix
        will be used.

    **kwargs : Passed to ``~scipy.ndimage.gaussian_filter``. The mode 
        default is "nearest" if none is given.
    
    Returns
    -------
    : ``~sunpy.map.Map``
        The altered Sunpy map object
    """
    gaussian_filter = {"mode":"nearest"} | kwargs
    gaussian_width_arcsec = FWHM_ARCSEC/(2*np.sqrt(2*np.log(2)))
    sigma = gaussian_width_arcsec.value/IMAGE_RES if sigma is None else sigma
    dd = ndimage.gaussian_filter(sunpy_map_obj.data, sigma, mode=gaussian_filter["mode"])
    return assign_plot_settings(sunpy.map.Map(dd, sunpy_map_obj.meta))

def deconvolve_with_array(sunpy_map_obj, psf_array, **kwargs):
    """Take a map and a point spread function and deconvolve using the 
    Richardson-Lucy method with a number of iterations. 

    [1] https://scikit-image.org/docs/stable/api/skimage.restoration.html#skimage.restoration.richardson_lucy

    Parameters
    ----------
    sunpy_map_obj : ``~sunpy.map.Map``
        The map for the filter to be applied to.

    psf_array : ``~numpy.ndarray``
        The PSF you want to use. 

    **kwargs : Passed to ``~skimage.restoration.richardson_lucy``. The
        Default for the clip input will be False and the number of 
        iterations (num_iter) will be 100.
    """
    deconvolve_filter = {"num_iter":100, "clip":False} | kwargs
    deconvolved = restoration.richardson_lucy(sunpy_map_obj.data, psf_array, **deconvolve_filter)
    return assign_plot_settings(sunpy.map.Map(deconvolved, sunpy_map_obj.meta))

def deconvolve_with_file(sunpy_map_obj, psf_file, oa_to_source=None, hor_to_source=None, **kwargs):
    """Take a map and a point spread function and deconvolve using the 
    Richardson-Lucy method with a number of iterations. 

    [1] https://scikit-image.org/docs/stable/api/skimage.restoration.html#skimage.restoration.richardson_lucy

    Parameters
    ----------
    sunpy_map_obj : ``~sunpy.map.Map``
        The map for the filter to be applied to.

    psf_file : `str`
        The file containing the NuSTAR PSF. This should be in the user's 
        CALDB install, in the location \".../caldb/data/nustar/fpm/bcf/psf/\"
        with the name format \"nu<FPM>2dpsfen1_20100101v001.fits\", where 
        <FPM> will be A or B.

    oa_to_source : float
        Angle subtended between the optical axis (OA), observer, and the 
        X-ray source in arcminutes (0<=oa_to_source<8.5 arcminutes), 
        i.e. radial distance to the source from the OA. Chooses the 
        correct PSF data to use. If None, then assumes source is on-axis 
        with an `oa_to_source` of 0.
        Default: None

    hor_to_source : float
        Angle subtended between horizontal through the optical axis (OA), 
        and the line through the X-ray source and OA in degrees.
        Clockwise is positive and anticlockwise is negative. Symmetric 
        reflected in the origin so -90<=hor_to_source<=90. If None, 
        assumes no rotation is needed with an `hor_to_source` of 0.
        Default: None

    **kwargs : Passed to ``~skimage.restoration.richardson_lucy``. The
        Default for the clip input will be False and the number of 
        iterations will be 10.
    """
    oa_to_source = 0 if oa_to_source is None else oa_to_source
    hor_to_source = 0 if hor_to_source is None else hor_to_source

    hdr_unit = _choose_os_angle_header_unit(oa_to_source)
    with fits.open(psf_file) as psf_hdulist:
        psf_res = psf_hdulist[hdr_unit].header['CDELT1'] # increment in degrees/pix
        psf_array = psf_hdulist[hdr_unit].data

    # check same res, at least in 1-D
    assert psf_res*3600 == sunpy_map_obj.meta['CDELT1'], "The resolution in the PSF and the current map are different."
    
    psf_array = ndimage.rotate(psf_array, hor_to_source, reshape=True)
    return deconvolve_with_array(sunpy_map_obj, psf_array, **kwargs)

def _choose_os_angle_header_unit(oa_to_source):
    """Use the OA-to-source angle to select the correct PSF."""
    # angles of 0 to 8.5 arcmin in 0.5 arcmin increments
    psf_oa_angles = np.arange(0,9,0.5) 
    # find the closest arcmin array
    index = np.argmin([abs(psfoaangles - oa_to_source) for psfoaangles in psf_oa_angles]) 
    return index+1