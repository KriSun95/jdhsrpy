from astropy.time import Time
import astropy.units as u
import matplotlib.colors

NUSTAR_EPOCH = Time("2010-01-01T00:00:00.000", format='isot',scale='utc') 
PIXELUNIT = "arcsec"
IMAGE_RES = 2.45810736
FWHM_ARCSEC = 18 << u.Unit(PIXELUNIT)
IMAGE_COLORMAP = "Spectral_r"
IMAGE_SCALE = matplotlib.colors.LogNorm()
