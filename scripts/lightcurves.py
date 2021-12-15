""" This set of functions is primarily used for photometery."""
from config import Configuration
from libraries.utils import Utils
from libraries.photometry import Photometry
from astropy.io import fits
from astropy.time import Time
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


class Lightcurves:

    @staticmethod
    def psf_phot(star):
        """ This function is the wrapper function for the aperture photometry.
        :parameter star - The star for the data to be reduced

        :return star_list - The data frame with the star information
        """

        # get the file list for the photometry
        files = Utils.get_file_list(Configuration.CLEAN_DIRECTORY, Configuration.FILE_EXTENSION)

        # make a data frame with the image information
        images = pd.DataFrame(columns=['filepath', 'jd', 'exptime'])
        lc = pd.DataFrame(columns=['jd', 'mag', 'er', 'x', 'y'])

        # loop through each image
        for idx, file in enumerate(files):

            # read in the image
            img = fits.getdata(Configuration.CLEAN_DIRECTORY + file)
            header = fits.getheader(Configuration.CLEAN_DIRECTORY + file)

            # do the aperture photometry for the stars
            star_list = Photometry.single_frame_psf_photometry(img,
                                                               Configuration.CLEAN_DIRECTORY + '/' + file.split('.')[0])

            if (idx % 50 == 0) & (idx > 0):
                Utils.log("Completed 50 images, working on the next 50. " +
                          str(len(files)-idx) + " files remain. ", "info")

        return

    @staticmethod
    def aperture_phot(star):
        """ This function is the wrapper function for the aperture photometry.
        :parameter star - The star for the data to be reduced

        :return star_list - The data frame with the star information
        """

        # get the file list for the photometry
        files = Utils.get_file_list(Configuration.CLEAN_DATE_DIRECTORY, Configuration.FILE_EXTENSION)

        # make a data frame with the image information
        images = pd.DataFrame(columns=['filepath', 'jd', 'exptime'])
        lc = pd.DataFrame(columns=['jd', 'mag', 'er', 'x', 'y'])

        # loop through each image
        for idx, file in enumerate(files):

            # read in the image
            img = fits.getdata(Configuration.CLEAN_DATE_DIRECTORY + file)
            header = fits.getheader(Configuration.CLEAN_DATE_DIRECTORY + file)

            # get the julian date of the image
            t = Time(header['DATE'])
            jd = t.jd + float(header['RELSEC']) / 3600. / 24.

            # update the image data frame
            images.loc[len(images.index)] = [Configuration.CLEAN_DATE_DIRECTORY + file, jd, header['EXPOSURE']]

            # get the stars on the image
            star_list = Photometry.single_frame_star_finder(img)

            # do the aperture photometry for the stars
            star_list = Photometry.single_frame_aper(star_list,
                                                     img,
                                                     Configuration.CLEAN_DATE_DIRECTORY + '/' + file.split('.')[0])

            lc.loc[len(lc.index)] = [jd, star_list.mag.min(), star_list.mag_err.min(),
                                     star_list[star_list.mag == star_list.mag.min()].x.min(),
                                     star_list[star_list.mag == star_list.mag.min()].y.min()]

            if (idx % 50 == 0) & (idx > 0):
                Utils.log("Completed 50 images, working on the next 50. " +
                          str(len(files)-idx) + " files remain. ", "info")

        return

    @staticmethod
    def mk_raw_lightcurves(flux_dir):
        """ This function will create the individual raw light curve files for each star in the specific star list

        :parameter flux_dir - A data frame with the master frame flux data

        :return nothing is returned, but each light curve is output
        """

        # get the file list of the differenced image flux information
        files = Utils.get_file_list(flux_dir, '.flux')

        # get the star list based on the x/y measurements
        star_list = Photometry.get_star_list(flux_dir, files)

        # combine the flux from the flux files, and write the raw light curves
        Photometry.combine_flux_files(flux_dir, files, star_list)

        return
