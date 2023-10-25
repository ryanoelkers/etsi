from libraries.utils import Utils
from scripts.master import Master
from scripts.clean import Clean
from libraries.photometry import Photometry
from config import Configuration

# do necessary prep work such as making output directories
Utils.create_directories(Configuration.DIRECTORIES)

if Configuration.PER_BAND == 'N':
    # generate the PSF if needed
    Utils.log("Generating co-adds and PSF for " + Configuration.STAR, "info")
    Master.generate_coadds(Configuration.COADD_IMAGE)

    # clean the images, if you want to
    if Configuration.CLEAN == 'Y':
        Utils.log("Cleaning images for " + Configuration.STAR, "info")
        Clean.clean_images(bias_subtract='Y')

    # generate the photometry if needed
    Utils.log("Starting Photometry for " + Configuration.STAR, "info")
    Photometry.psf_and_aperture_photometry()

else:
    Photometry.single_band_photometry(Configuration.COADD_DIRECTORY, Configuration.BAND)

Utils.log("All done! See ya later, alligator.", "info")
