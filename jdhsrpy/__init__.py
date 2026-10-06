import astropy.units as u
import matplotlib.colors
from astropy.time import Time

TEST_DATA_LOCATION = "http://foxsi.space.umn.edu/data/tmp/jdhsr"

NUSTAR_EPOCH = Time("2010-01-01T00:00:00.000", format="isot", scale="utc")
PIXELUNIT = "arcsec"
IMAGE_RES = 2.45810736
FWHM_ARCSEC = 18 << u.Unit(PIXELUNIT)
IMAGE_COLORMAP = "Spectral_r"
IMAGE_SCALE = matplotlib.colors.LogNorm()
