import astropy.units as u
import matplotlib.colors
from astropy.time import Time

TEST_DATA_LOCATION = "https://foxsi.space.umn.edu/JDHSR"

NUSTAR_EPOCH_ISOT_UTC = Time("2010-01-01T00:00:00.000", format="isot", scale="utc")
NUSTAR_EPOCH_MJD = ((55197 << u.day) - (5 << u.second)).to(u.day)
NUSTAR_EPOCH = Time(NUSTAR_EPOCH_MJD, format="mjd").utc
NUSTAR_EPOCH.format = "isot"
PIXELUNIT = "arcsec"
IMAGE_RES = 2.45810736 << u.arcsec / u.pixel
FWHM_ARCSEC = 18 << u.Unit(PIXELUNIT)
IMAGE_COLORMAP = "Spectral_r"
IMAGE_SCALE = matplotlib.colors.LogNorm()
RAW_PIX_DIFF = 12.646031412008877 << u.arcsec / u.pixel
CHIP_GAP = 15.77588557692691 << u.arcsec
PIX_GRID = (32, 32) << u.pixel
