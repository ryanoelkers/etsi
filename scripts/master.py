""" This set of functions is primarily used for photometery."""
from config import Configuration
from libraries.utils import Utils
from libraries.photometry import Photometry
import numpy as np
import warnings
import os
import astroalign as aa
from astropy.io import fits
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)

class Master:
    @staticmethod
    def generate_master_files(image_directory, combine_type='median'):
        """ This script will generate the master frame and generate the PSF to use for photometry.

        :parameter image_directory - The directory with the cleaned images for the data reduction
        :parameter combine_type - The type of combination for the images, if not median.

        :return - Nothing is returned, but the master frame is generated and the PSf created
        """

        # generate, or return, the master frame
        Master.mk_master(image_directory)

        # generate the PSF required for the image
        epsf = Photometry.generate_psf(Configuration.MASTER_DIRECTORY)

        return epsf

    @staticmethod
    def mk_master(image_directory, combine_type='median'):
        """ This function will make the master frame using the provided image list.
        :parameter image_directory - a directory where the images reside for combination
        :parameter combine_type - Either median or mean depending on how you want to combine the files

        :return - The bias frame is returned and written to the master directory
        """

        if os.path.isfile(Configuration.MASTER_DIRECTORY + Configuration.STAR + '_master.fits') == 0:

            if combine_type == 'median':
                # get the image list
                image_list = Utils.get_file_list(image_directory, Configuration.FILE_EXTENSION)

                # determine the number of loops we need to move through for each image
                nfiles = len(image_list)
                nbulk = 100

                # get the integer and remainder for the combination
                full_bulk = nfiles // nbulk
                part_bulk = nfiles % nbulk
                if part_bulk > 0:
                    hold_bulk = full_bulk + 1
                else:
                    hold_bulk = full_bulk

                # here is the 'holder'
                hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXS_X, Configuration.AXS_Y))

                # update the log
                Utils.log("Generating a master frame from multiple files in bulks of " + str(nbulk) +
                          " images. There are " + str(nfiles) + " images to combine, which means there should be " +
                          str(hold_bulk) + " mini-files to combine.", "info")
                ref_img = fits.getdata(image_directory + image_list[0])
                if np.ndim(ref_img) > 2:
                   ref_img = ref_img[0]

                for kk in range(0, hold_bulk):

                    # loop through the images in sets of nbulk
                    if kk < full_bulk:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(nbulk, Configuration.AXS_X, Configuration.AXS_Y))

                        # generate the max index
                        mx_index = nbulk
                    else:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(part_bulk, Configuration.AXS_X, Configuration.AXS_Y))

                        # generate the max index
                        mx_index = part_bulk

                    # make the starting index
                    loop_start = kk * nbulk
                    idx_cnt = 0

                    Utils.log("Making mini file " + str(kk) + ".", "info")

                    # now loop through the images
                    for jj in range(loop_start, mx_index + loop_start):
                        # read in the image directly into the block_hold
                        img = fits.getdata(image_directory + image_list[jj])

                        if np.ndim(img) > 2:
                            img = img[0]

                        try:
                            master_tmp, footprint = aa.register(ref_img.byteswap().newbyteorder(),
                                                                img.byteswap().newbyteorder(),
                                                                min_area=9)
                        except aa.MaxIterError:
                            master_tmp = img
                        except ValueError:
                            master_tmp = img
                        block_hold[idx_cnt] = master_tmp

                        # increase the iteration
                        idx_cnt += 1

                    # median the data into a single file
                    hold_data[kk] = np.median(block_hold, axis=0)

                # median the mini-images into one large image
                master = np.median(hold_data, axis=0)

                # pull the header information from the first file of the set
                master_header = fits.getheader(image_directory + image_list[0])

                master_header['MAST_COMB'] = 'median'
                master_header['NUM_MAST'] = nfiles

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + Configuration.STAR + '_master.fits',
                             master, master_header, overwrite=True)

        else:
            Utils.log('Master frame found. Using legacy file, please delete if not wanted!', 'info')

        return