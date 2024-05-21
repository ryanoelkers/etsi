""" This set of functions is primarily used for photometery."""
from config import Configuration
from libraries.utils import Utils
from libraries.photometry import Photometry
from libraries.preprocessing import Preprocessing
from scripts.clean import Clean
import numpy as np
import warnings
import os
import astroalign as aa
from astropy.io import fits
from astropy import time
from dateutil import parser
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)


class Master:

    @staticmethod
    def generate_coadds(coadd='Y'):
        """ This script will generate co-adds with the intention that photometry should be run on high SNR frames

        :parameter coadd = Y/N if you want to generated co-adds

        :return - Nothing is returned, but the co-adds are generated
        """

        if coadd == 'Y':
            # generate the number of images to combine based on the bin time and the exposure time of the observations
            Master.mk_coadds(Configuration.RAW_DIRECTORY,
                             Configuration.BIN_NUM,
                             Configuration.BIN_TIME,
                             combine_type='time')
        else:
            # generate the PSF required for the image
            Utils.log("Skipping co-adds.", 'info')

        return

    @staticmethod
    def generate_master_files(image_directory, num_to_combine=1000):
        """ This script will generate the master frame and generate the PSF to use for photometry.

        :parameter image_directory - The directory with the cleaned images for the data reduction
        :parameter num_to_combine - The number of images to combine for the master frame, if not 1000

        :return - the PSF is returned
        """

        # generate, or return, the master frame
        if Configuration.MAKE_MASTER == 'Y':
            Master.mk_master(image_directory, num_to_combine)
        else:
            Utils.log('Skipping generation of Master Frame. Please change Y/N if code breaks!', 'info')

        # generate the PSF required for the image
        epsf = Photometry.generate_psf(Configuration.MASTER_DIRECTORY, find_stars='N', plot_psf='N')

        return epsf

    @staticmethod
    def mk_master(image_directory, num_to_combine, combine_type='median'):
        """ This function will make the master frame using the provided image list.
        :parameter image_directory - a directory where the images reside for combination
        :parameter combine_type - Either median or mean depending on how you want to combine the files
        :parameter num_to_combine - The number of images to comnbine
        :return - The bias frame is returned and written to the master directory
        """

        if Configuration.SKY_SUBTRACT == 'global':
            mst_nme = '_master_no_bkg.fits'
        else:
            mst_nme = '_master.fits'

        if os.path.isfile(Configuration.MASTER_DIRECTORY +
                          Configuration.STAR + '_' + Configuration.BEAM_TYPE + mst_nme) == 0:

            if combine_type == 'median':
                # get the image list
                image_list = Utils.get_file_list(image_directory, '.fits')

                # determine the number of loops we need to move through for each image
                if len(image_list) < 100:
                    nfiles = len(image_list)
                    nbulk = len(image_list)
                else:
                    nfiles = num_to_combine
                    nbulk = 100

                # get the integer and remainder for the combination
                full_bulk = nfiles // nbulk
                part_bulk = nfiles % nbulk
                if part_bulk > 0:
                    hold_bulk = full_bulk + 1
                else:
                    hold_bulk = full_bulk

                # here is the 'holder'
                hold_data = np.ndarray(shape=(hold_bulk, Configuration.AXIS_X, Configuration.AXIS_Y))

                # update the log
                Utils.log("Generating a master frame from multiple files in bulks of " + str(nbulk) +
                          " images. There are " + str(nfiles) + " images to combine, which means there should be " +
                          str(hold_bulk) + " mini-files to combine.", "info")

                # clean the image
                ref_img, ref_header = fits.getdata(image_directory + image_list[0], 0, header=True)

                for kk in range(0, hold_bulk):

                    # loop through the images in sets of nbulk
                    if kk < full_bulk:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(nbulk, Configuration.AXIS_X, Configuration.AXIS_Y))

                        # generate the max index
                        mx_index = nbulk
                    else:
                        # generate the image holder
                        block_hold = np.ndarray(shape=(part_bulk, Configuration.AXIS_X, Configuration.AXIS_Y))

                        # generate the max index
                        mx_index = part_bulk

                    # make the starting index
                    loop_start = kk * nbulk
                    idx_cnt = 0

                    Utils.log("Making mini file " + str(kk) + ".", "info")

                    # now loop through the images
                    for jj in range(loop_start, mx_index + loop_start):
                        # read in the image directly into the block_hold

                        # clean the image
                        img, header = fits.getdata(image_directory + image_list[jj], 0, header=True)

                        if np.ndim(img) > 2:
                            img = img[0]

                        try:
                            master_tmp, footprint = aa.register(ref_img.byteswap().newbyteorder(),
                                                                img.byteswap().newbyteorder(), max_control_points=50,
                                                                detection_sigma=10, min_area=50)
                        except aa.MaxIterError:
                            master_tmp = img
                        except ValueError:
                            master_tmp = img
                        if Configuration.SKY_SUBTRACT == 'global':
                            master_tmp, master_tmp_header = Preprocessing.sky_subtract(master_tmp,
                                                                                       header,
                                                                                       sky_write='N')
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
                fits.writeto(Configuration.MASTER_DIRECTORY +
                             Configuration.STAR + '_' + Configuration.BEAM_TYPE + mst_nme,
                             master, master_header, overwrite=True)
        else:
            Utils.log('Master frame found. Using legacy file, please delete if not wanted!', 'info')

        return

    @staticmethod
    def mk_coadds(image_directory, num_to_combine, time_to_combine, combine_type='time'):

        """ This function will make the master frame using the provided image list.
        :parameter image_directory - a directory where the images reside for combination
        :parameter num_to_combine - The number of images to combine
        :parameter time_to_combine - The approximate exposure time (in seconds) to combine the data
        :parameter combine_type - either time or number to indicate how to co-add the images

        :return - The bias frame is returned and written to the master directory
        """
        # get the image list
        image_list = Utils.get_file_list(image_directory, Configuration.FILE_EXTENSION)

        # determine the number of loops we need to move through for each image
        nfiles = len(image_list)

        # update the log
        Utils.log("Generating co-add frames from multiple files. There are " + str(nfiles) + " images to " +
                  "co-add.", "info")

        coadd_flag = 0
        img_idx = 0

        for kk in range(0, nfiles):

            if coadd_flag == 0:
                # here is the 'holder'
                hold_data = np.zeros(shape=(Configuration.AXIS_X, Configuration.AXIS_Y))
                idx_cnt = 0
                img_exp = 0

            # loop through the images in sets of nbulk
            if kk % 100 == 0:
                Utils.log("Making co-add file " + str(img_idx) + ". There are " +
                          str(nfiles - kk) + " left for co-adds.", "info")

            # clean the image
            try:
                img, header, bd_flag = Clean.calibrate_img(image_list[kk], image_directory,
                                                           bias_subtract=Configuration.BIAS_SUBTRACT)
            except:
                Utils.log("Bad image is: " + image_list[kk], "info")
                continue

            if np.ndim(img) > 2:
                img = img[0]

            # get the time stamps for the co-addition
            if coadd_flag == 0:

                # get the time start point for the co-add image
                jd_st = time.Time(parser.parse(header['DATE-OBS'])).jd
                # get the airmass start point for the co-add image
                airmass_st = header['AIRMASS']
                if Configuration.FILE_EXTENSION == '.fits':
                    focus_st = header['FOCUS']
                    press_st = header['HIERARCH PRESSURE_MB']
                    temp_st = header['HIERARCH TEMPERATURE_C']
                    humid_st = header['HIERARCH HUMIDITY_PERCENT']
                    # dome_st = header['HIERARCH DOME_HOUR']
                    zenith_st = header['HIERARCH ZD_DEGREE']
                    azimuth_st = header['HIERARCH AZ_DEGREE']
                # some data did not have the exposure time properly in the header, if not, make an approximate exp time
                try:
                    if (Configuration.BEAM_TYPE == 'reflection') & (Configuration.FILE_EXTENSION == '.fits'):
                        img_exp = img_exp + ((header['HIERARCH EXPOSURE TIME'] / 1000) * int(header['CO-ADD']))
                    else:
                        img_exp = img_exp + (np.around(header['HIERARCH EXPOSURE TIME'], decimals=2) * header['CO-ADD'])
                except KeyError:
                    img_exp = img_exp + (Configuration.EXPOSURE_TIME * Configuration.EXPOSURE_COADD)

                # co-add the image & increase the image count
                hold_data = hold_data + img
                idx_cnt = idx_cnt + 1

            if coadd_flag == 1:
                # update the "final" header time and airmass
                jd_ed = time.Time(parser.parse(header['DATE-OBS'])).jd
                airmass_ed = header['AIRMASS']
                if Configuration.FILE_EXTENSION == '.fits':
                    focus_ed = header['FOCUS']
                    press_ed = header['HIERARCH PRESSURE_MB']
                    temp_ed = header['HIERARCH TEMPERATURE_C']
                    humid_ed = header['HIERARCH HUMIDITY_PERCENT']
                    # dome_ed = header['HIERARCH DOME_HOUR']
                    zenith_ed = header['HIERARCH ZD_DEGREE']
                    azimuth_ed = header['HIERARCH AZ_DEGREE']
                # some data did not have the exposure time properly in the header, if not, make an approximate exp time
                try:
                    if (Configuration.BEAM_TYPE == 'reflection') & (Configuration.FILE_EXTENSION == '.fits'):
                        img_exp = img_exp + ((header['HIERARCH EXPOSURE TIME'] / 1000) * int(header['CO-ADD']))
                    else:
                        img_exp = img_exp + (np.around(header['HIERARCH EXPOSURE TIME'], decimals=2) * header['CO-ADD'])
                except KeyError:
                    img_exp = img_exp + (Configuration.EXPOSURE_TIME * Configuration.EXPOSURE_COADD)

                # co-add the image & increase the image count
                hold_data = hold_data + img
                idx_cnt = idx_cnt + 1

            if coadd_flag == 1:
                time_offset_check = (jd_ed - jd_st) * 24. * 60. * 60
            else:
                time_offset_check = 0.

            # check if co-add size has been reached (either from time or number)
            if ((combine_type == 'time') & (img_exp >= time_to_combine)) | \
                    ((combine_type == 'number') & (idx_cnt >= num_to_combine)) | \
                    (kk >= nfiles-1) | (time_offset_check > 5 * time_to_combine):
                # pull the header information from the first file of the set
                hold_header = fits.PrimaryHDU().header  # fits.getheader(image_directory + image_list[0])

                # update the header with the co-add information
                hold_header['COADD_NUM'] = idx_cnt
                hold_header['COADD_ST'] = np.around(jd_st, decimals=6)
                hold_header['COADD_ED'] = np.around(jd_ed, decimals=6)
                hold_header['EXP_TIME'] = img_exp
                hold_header['AIRMASSS'] = airmass_st
                hold_header['AIRMASSE'] = airmass_ed
                if Configuration.FILE_EXTENSION == '.fits':
                    hold_header['HUMIDS'] = humid_st
                    hold_header['HUMIDE'] = humid_ed
                    hold_header['FOCUSS'] = focus_st
                    hold_header['FOCUSE'] = focus_ed
                    hold_header['PRESSS'] = press_st
                    hold_header['PRESSE'] = press_ed
                    hold_header['TEMPS'] = temp_st
                    hold_header['TEMPE'] = temp_ed
                    # hold_header['DOMES'] = dome_st
                    # hold_header['DOMEE'] = dome_ed
                    hold_header['ZENITHS'] = zenith_st
                    hold_header['ZENITHE'] = zenith_ed
                    hold_header['AZIMS'] = azimuth_st
                    hold_header['AZIME'] = azimuth_ed

                # write the image out to the master directory
                if img_idx < 10:
                    file_numval = '000' + str(img_idx)
                elif (img_idx >= 10) & (img_idx < 100):
                    file_numval = '00' + str(img_idx)
                elif (img_idx >= 100) & (img_idx < 1000):
                    file_numval = '0' + str(img_idx)
                else:
                    file_numval = str(img_idx)

                # write the new fits file to the co-add directory
                if combine_type == 'number':
                    hold_data = hold_data / Configuration.BIN_NUM
                fits.writeto(Configuration.COADD_DIRECTORY +
                             Configuration.STAR + '_' + Configuration.BEAM_TYPE + '_' + file_numval + '.fits',
                             hold_data, hold_header, overwrite=True)

                Utils.log('Co-add #' + str(img_idx) + ' generated.', "info")

                # reset counters
                img_idx = img_idx + 1
                idx_cnt = 0
                img_exp = 0
                coadd_flag = 0
            else:
                coadd_flag = 1
                continue
        return
