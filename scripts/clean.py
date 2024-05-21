""" This script will clean the FFI images, it has been significantly reduced to only
show specific steps. The majority of the functional library can be found in the libraries
directory in utils, fits, photometry, image, and preprocessing libraries."""
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
from config import Configuration
import os
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import datetime
import pytz
import time
import numpy as np
from skimage import io
from PIL import Image
import logging
logging.getLogger("PIL").setLevel(logging.WARNING)


class Clean:

    @staticmethod
    def clean_images(sky_subtract="N", bias_subtract="N", flat_divide='N', alignment='N', dark_subtract="N"):
        """ This is the main function script to clean multiple images, alternatively clean_img can be used to clean
        a single image.
        :parameter sky_subtract - Y/N if you want to subtract the sky from the images (default = Y)
        :parameter bias_subtract - Y/N if you want to subtract the bias from the images (default = N)
        :parameter flat_divide - Y/N if you want to flatten the images (default = N)
        :parameter alignment - Y/N if you want to align the images (default = N)
        :parameter dark_subtract -Y/N if you want to subtract a dark frame from the science images (default = N)

        return no value is returned, the values images from in_path are cleaned and deposited in out_path
        """
        st = time.time()  # clock started

        # get the file list
        Utils.log("Getting file list...", "info")
        if Configuration.COADDED == 'Y':
            img_dir = Configuration.COADD_DIRECTORY
        else:
            img_dir = Configuration.RAW_DIRECTORY

        files = Utils.get_file_list(img_dir, Configuration.FILE_EXTENSION)

        # break if there are no files
        if len(files) == 0:
            Utils.log("No .fits files found in " + img_dir + ". Breaking...", "debug")
            return()
        
        Utils.log("Starting to clean " + str(len(files)) + " images.", "info")
        for idx, file in enumerate(files):

            # make a new name for the file based on which actions are taken
            file_name = Preprocessing.mk_nme(file, 'N', sky_subtract, bias_subtract,
                                             flat_divide, alignment, dark_subtract)

            # only create the files that don't exist
            if os.path.isfile(Configuration.CLEAN_DIRECTORY + file_name) == 1:
                Utils.log("Image " + Configuration.CLEAN_DIRECTORY + file_name +
                          " already exists. Skipping for now...", "info")

            # if the image does not exist then clean
            if os.path.isfile(Configuration.CLEAN_DIRECTORY + file_name) == 0:

                # clean the image
                clean_img, header, bd_flag = Clean.calibrate_img(file, img_dir, sky_subtract,
                                                                 bias_subtract, flat_divide, alignment, dark_subtract)

                # write out the file
                if bd_flag == 0:
                    fits.writeto(Configuration.CLEAN_DIRECTORY + file_name,
                                 clean_img, header, overwrite=True)

                    # print an update to the cleaning process
                    # Utils.log("Cleaned image written as " +
                    #          Configuration.CLEAN_DIRECTORY + file_name + ".", "info")
                else:
                    Utils.log(file_name + " is a bad image. Not written.", "info")

            Utils.log(str(len(files) - idx - 1) + " images remain to be cleaned.",  "info")

        fn = time.time()  # clock stopped
        Utils.log("Imaging cleaning complete in " + str(np.around((fn - st), decimals=2)) + "s.", "info")

    @staticmethod
    def calibrate_img(file, ref_path, sky_subtract='N', bias_subtract='N', flat_divide='N', dark_subtract="N"):
        """ This function is the primary script to clean the image, various other functions found in this class
        can be found in the various libraries imported.

        :parameter  file - The file name of the image you would like to clean
        :parameter ref_path - The path to the reference frame
        :parameter sky_subtract - Y/N if you want to do the global sky subtraction
        :parameter bias_subtract - Y/N if you want to remove a bias frame (default = N)
        :parameter flat_divide - Y/N if you want to flatten the image (default = N)
        :parameter dark_subtract -Y/N if you want to subtract the dark frame (default = N)
        """

        # Utils.log("Now cleaning " + file + ".", "info")

        # read in the image
        if Configuration.FILE_EXTENSION == '.fits':
            img, header = fits.getdata(ref_path + file, 0, header=True)
        else:
            img = io.imread(ref_path + file)

            # generate the header with the time stamp, exposure, and airmass
            header = fits.PrimaryHDU().header

            # time stamp
            dte = datetime.datetime.fromtimestamp(os.path.getmtime(ref_path + file))
            local = pytz.timezone("America/Chicago")
            local_dte = local.localize(dte, is_dst=Configuration.DST)
            utc_dte = local_dte.astimezone(pytz.utc)
            # if Configuration.MONTH != 'July':
            header['DATE-OBS'] = str(utc_dte)
            # else:
            #     header['DATE-OBS'] = str(dte)

            # airmass
            observatory = EarthLocation(lat=Configuration.OBS_LAT * u.deg,
                                        lon=Configuration.OBS_LON * u.deg,
                                        height=Configuration.OBS_HGT * u.m)
            obj_coord = SkyCoord(ra=Configuration.RA * u.deg,
                                 dec=Configuration.DEC * u.deg,
                                 frame='icrs')
            obj_altaz = obj_coord.transform_to(AltAz(obstime=utc_dte, location=observatory))
            header['AIRMASS'] = np.around(obj_altaz.secz.value, decimals=6)

            # exposure
            with Image.open(ref_path + file) as img_head:

                header['HIERARCH EXPOSURE TIME'] = int(''.join(filter(lambda x: x.isdigit(),
                                                                      img_head.tag_v2[270].
                                                                      split("\n")[5].split("=")[1]))) / 1000
            header['CO-ADD'] = 1
        if np.ndim(img) > 2:
            img = img[0]

        # bias subtract if necessary
        if bias_subtract == 'Y':
            bias = Preprocessing.mk_bias(Configuration.BIAS_DIRECTORY, dark='N', combine_type=Configuration.BIAS_TYPE)
            img, header = Preprocessing.bias_subtract(img, header, dark='N')

        # dark subtract if necessary
        if dark_subtract == 'Y':
            dark = Preprocessing.mk_bias(Configuration.DARKS_DIRECTORY, dark='Y', combine_type='median')
            img, header = Preprocessing.bias_subtract(img, header, dark='Y')

        # flat divide if necessary
        if flat_divide == 'Y':
            flat = Preprocessing.mk_flat(Configuration.FLATS_DIRECTORY)
            img, header = Preprocessing.flat_divide(img, header)

        if sky_subtract == 'Y':
            img, header = Preprocessing.sky_subtract(img, header, sky_write='Y')

        # Utils.log("Cleaning finished.", "info")
        bd_flag = 0

        return img, header, bd_flag
