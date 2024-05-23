""" This set of functions is primarily used for co-adding ETSI images."""
from config import Configuration
from libraries.utils import Utils
import astroalign as aa
import numpy as np
import warnings
from astropy.io import fits
import os
from skimage import io
from PIL import Image
from astropy import time
from dateutil import parser
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import datetime
import pytz
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)


class Calibration:

    @staticmethod
    def calibrate_to_start(img0, img, og_list, pv_list):
        """ This function will calibrate the img file to the img0 file

        :parameter img0 - The image to calibrate positions to
        :parameter img - The image you want to calibrate
        :parameter og_list - The list of the OG pixel positions
        :parameter pv_list - The previous list of positions

        :return star_list, flag - The updated star positions based on the first frame's data and 0/1 if it worked
        """

        # set up the coordinate list for the new image
        src = list()
        star_list = og_list.copy().reset_index(drop=True)

        # pull out each star pair
        for idy, row in star_list.iterrows():
            ppair = (row.x, row.y)
            src.append(ppair)

        src = np.array(src)
        star_list = star_list.copy().reset_index(drop=True)

        # set up the image transformation
        try:
            # get the transformation offset between the frames; these parameters may need to be updated!
            img_transf, (i0_list, i1_list) = aa.find_transform(img0 - np.nanmedian(img0),
                                                               img - np.nanmedian(img),
                                                               max_control_points=3,
                                                               detection_sigma=10)
            # get the new positions
            img_calc = aa.matrix_transform(src, img_transf.params)

            # update the star list with the positions
            star_list['x'] = img_calc[:, 0]
            star_list['y'] = img_calc[:, 1]

            # return the star list
            return star_list, 0

        # if the alignment fails, return the original positions
        except aa.MaxIterError:
            return pv_list, 1

        # if the alignment fails, return the original positions
        except ValueError:
            return pv_list, 1

    @staticmethod
    def generate_coadds(coadd='Y'):
        """ This script will generate co-adds with the intention that photometry should be run on high SNR frames

            :parameter coadd = Y/N if you want to generate co-adds

            :return - Nothing is returned, but the co-adds are generated
        """

        if coadd == 'Y':
            # if you want to co-add frame, then co-add the frames
            Calibration.mk_coadds(Configuration.RAW_DIRECTORY, Configuration.BIN_TIME)
        else:
            # if you don't want to co-add the frames, then do nothing
            Utils.log("You marked 'N' in the configuration file. No co-adds will be generated.", 'info')

        return

    @staticmethod
    def read_raw_etsi_data(image, directory):
        """This function will determine the ETSI data type and read in the image appropriately.

        :parameter image - The string with the image name
        :parameter directory - The string with the directory where the image is located

        :return img, header - The image and its header are returned.
        """

        # read in the image
        if Configuration.FILE_EXTENSION == '.fits':
            # if the image is a .fits file, we can just use astropy
            img, header = fits.getdata(directory + image, 0, header=True)
        else:
            # if the image is not a .fits file, we need to be a bit more involved
            img = io.imread(directory + image)

            # generate the header to hold the relevant image information
            header = fits.PrimaryHDU().header

            # the time stamp needs to be generated from the file creation time
            dte = datetime.datetime.fromtimestamp(os.path.getmtime(directory + image))

            # the data was all taken in Central Time, convert to GMT
            local = pytz.timezone("America/Chicago")
            local_dte = local.localize(dte, is_dst=Configuration.DST)

            utc_dte = local_dte.astimezone(pytz.utc)
            header['DATE-OBS'] = str(utc_dte)

            # calculate the expected airmass at McDonald
            observatory = EarthLocation(lat=Configuration.OBS_LAT * u.deg,
                                        lon=Configuration.OBS_LON * u.deg,
                                        height=Configuration.OBS_HGT * u.m)
            obj_coord = SkyCoord(ra=Configuration.RA * u.deg,
                                 dec=Configuration.DEC * u.deg,
                                 frame='icrs')
            obj_altaz = obj_coord.transform_to(AltAz(obstime=utc_dte, location=observatory))
            header['AIRMASS'] = np.around(obj_altaz.secz.value, decimals=6)

            # calculate the expected exposure time
            with Image.open(directory + image) as img_head:

                header['HIERARCH EXPOSURE TIME'] = int(''.join(filter(lambda x: x.isdigit(),
                                                                      img_head.tag_v2[270].
                                                                      split("\n")[5].split("=")[1]))) / 1000

            # this image was not co-added, so set the co-add flag to be 1
            header['CO-ADD'] = 1

        # make sure the image is 2-dimensional
        if np.ndim(img) > 2:
            img = img[0]

        return img, header

    @staticmethod
    def mk_coadds(image_directory, time_to_combine):
        """ This function will make co-adds based on the given time range to combine the files.

            :parameter image_directory - The directory where the images reside for co-addition.
            :parameter time_to_combine - The approximate exposure time (in seconds) you want to combine the data to.

            :return - Nothing is returned, but the co-added frame is written to the co-add directory in teh config file.
        """
        # get the image list
        image_list = Utils.get_file_list(image_directory, Configuration.FILE_EXTENSION)

        # get the number of files in the give directory
        nfiles = len(image_list)

        # provide information about the number of co-adds to be made
        Utils.log("Generating co-add frames from multiple files. There are " + str(nfiles) + " images to " +
                  "co-add.", "info")

        # initialize flags for dumping later
        coadd_flag = 0
        img_idx = 0

        # begin to iterate through all files
        for kk, image in enumerate(image_list):

            if coadd_flag == 0:
                # generate an empty array to use as the image
                hold_data = np.zeros(shape=(Configuration.AXIS_X, Configuration.AXIS_Y))

                # reset the image counter and expsoure counter
                idx_cnt = 0
                img_exp = 0

            # loop through the images in sets nbulk images to avoid overloading memory
            if kk % 100 == 0:
                Utils.log("Making co-add file " + str(img_idx) + ". There are " +
                          str(nfiles - kk) + " left for co-adds.", "info")

            # Some images are corrupted, if this is the case -- skip the image
            try:
                # read in the image if you can
                img, header = Calibration.read_raw_etsi_data(image, Configuration.RAW_DIRECTORY)
            except:
                # otherwise print to log file and skip
                Utils.log("Bad image is: " + image, "info")
                continue

            # if this is the first file in the sequence, initialize all of the header information
            if coadd_flag == 0:

                # get the initial values for the header
                jd_st = time.Time(parser.parse(header['DATE-OBS'])).jd
                airmass_st = header['AIRMASS']

            else:
                # update the "final" header time and airmass, it may not actually be final!
                jd_ed = time.Time(parser.parse(header['DATE-OBS'])).jd
                airmass_ed = header['AIRMASS']

            # some headers provide the wrong expsoure time, in these cases fix it!
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

            # sometimes there are gaps in the data due to the instrument or weather, if so make the coadd now!
            if coadd_flag == 1:
                # convert the difference in time to seconds
                time_offset_check = (jd_ed - jd_st) * 24. * 60. * 60
            else:
                # this is the first file so there is no offset
                time_offset_check = 0.

            # check if co-add size has been reached (either from time or pausing exposures)
            if (img_exp >= time_to_combine) | (kk >= nfiles-1) | (time_offset_check > 5 * time_to_combine):

                # make the co-added header
                hold_header = fits.PrimaryHDU().header

                # update the header with the co-add information
                hold_header['COADD_NUM'] = idx_cnt
                hold_header['COADD_ST'] = np.around(jd_st, decimals=6)
                hold_header['COADD_ED'] = np.around(jd_ed, decimals=6)
                hold_header['EXP_TIME'] = img_exp
                hold_header['AIRMASSS'] = airmass_st
                hold_header['AIRMASSE'] = airmass_ed

                # create the co-add ID number for the naming convention
                if img_idx < 10:
                    file_numval = '000' + str(img_idx)
                elif (img_idx >= 10) & (img_idx < 100):
                    file_numval = '00' + str(img_idx)
                elif (img_idx >= 100) & (img_idx < 1000):
                    file_numval = '0' + str(img_idx)
                else:
                    file_numval = str(img_idx)

                # write the new fits file to the co-add directory
                fits.writeto(Configuration.COADD_BEAM_DIRECTORY +
                             Configuration.STAR + '_' + Configuration.BEAM_TYPE + '_' + file_numval + '.fits',
                             hold_data, hold_header, overwrite=True)

                Utils.log('Co-add #' + str(img_idx) + ' generated.', "info")

                # reset counters
                img_idx = img_idx + 1
                idx_cnt = 0
                img_exp = 0
                coadd_flag = 0
            else:
                # if you are not at the time limit, then set coadd_flag to 1
                coadd_flag = 1

                continue
        return
