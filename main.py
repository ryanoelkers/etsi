from libraries.utils import Utils
from scripts.calibration import Calibration
from scripts.photometry import Photometry
from config import Configuration

# do necessary prep work such as making output directories
Utils.create_directories(Configuration.DIRECTORIES)

# generate the PSF if needed
Utils.log("Generating co-adds for " + Configuration.STAR, "info")
Calibration.generate_coadds(Configuration.COADD_IMAGE)

# generate the photometry if needed
Utils.log("Starting Photometry for " + Configuration.STAR, "info")
Photometry.run_photometry(Configuration.COADD_BEAM_DIRECTORY, Configuration.WVE_STR)

Utils.log("All done! See ya later, alligator.", "info")
