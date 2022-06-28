from libraries.utils import Utils
from scripts.lightcurves import Lightcurves
from scripts.master import Master
from scripts.clean import Clean
from libraries.photometry import Photometry
from config import Configuration

# do necessary prep work such as making output directories
Utils.create_directories(Configuration.DIRECTORIES)

# do the necessary preprocessing of the images
if Configuration.CLEAN_SKIP == 'N':
    Clean.clean_images(bias_subtract=Configuration.BIAS_SUBTRACT, sky_subtract=Configuration.SKY_SUBTRACT,
                       dark_subtract=Configuration.DARK_SUBTRACT, flat_divide=Configuration.FLAT_DIVIDE,
                       alignment=Configuration.ALIGNMENT)
else:
    Utils.log("Skipping image cleaning.", "info")

# run the photometry on the images
if Configuration.PHOTOMETRY_SKIP == 'N':

    # start the photometry
    epsf = Master.generate_master_files(Configuration.CLEAN_DIRECTORY)
    Photometry.psf_and_aperture_photometry(epsf)

else:
    Utils.log("Skipping photometry generation.", "info")

if Configuration.LIGHTCURVE_SKIP == 'N':

    # detrend the photometry
    Lightcurves.detrend_lightcurves(Configuration.LIGHTCURVE_DIRECTORY + "raw\\")
else:
    Utils.log("Skipping ligth curve detrending.", "info")

Utils.log("All done! See ya later, alligator.", "info")
