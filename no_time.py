from libraries.utils import Utils
from scripts.master import Master
from libraries.photometry import Photometry
from config import Configuration

# do necessary prep work such as making output directories
Utils.create_directories(Configuration.DIRECTORIES)

# generate the PSF if needed
Utils.log("Generating PSF for " + Configuration.STAR, "info")
epsf = Master.generate_master_files(Configuration.RAW_DIRECTORY, Configuration.NUM_MASTER_FILES)

# generate the photometry if needed
Utils.log("Starting Photometry for " + Configuration.STAR, "info")
Photometry.psf_and_aperture_photometry(epsf)

Utils.log("All done! See ya later, alligator.", "info")
