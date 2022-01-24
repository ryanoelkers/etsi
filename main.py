from libraries.utils import Utils
from scripts.lightcurves import Lightcurves
from scripts.clean import Clean
from scripts.calibration import Calibration
from config import Configuration

# do necessary prep work such as making output directories
Utils.create_directories(Configuration.DIRECTORIES)

# do the necessary preprocessing of the images
if Configuration.CLEAN_SKIP == 'N':
    Clean.clean_images(dark_subtract='Y', sky_subtract='Y')
else:
    Utils.log("Skipping image cleaning.", "info")

# run the photometry on the images
if Configuration.PHOTOMETRY_SKIP == 'N':
    if Configuration.PHOTOMETRY == 'PSF':
        Lightcurves.psf_phot(Configuration.STAR)
    else:
        Lightcurves.aperture_phot(Configuration.STAR)
else:
    Utils.log("Skipping photometry generation.", "info")

if Configuration.LIGHTCURVE_SKIP == 'N':
    Lightcurves.mk_raw_lightcurves(Configuration.CLEAN_DIRECTORY)
else:
    Utils.log("Skipping light curve generation.", "info")

if Configuration.COLOR_SKIP == 'N':
    Lightcurves.color_comparison(Configuration.LIGHTCURVE_DIRECTORY)
else:
    Utils.log("Skipping color comparison.", "info")
# output the light curves
Utils.log("All done! See ya later, alligator.", "info")
