""" This set of functions is primarily used for photometery."""
import astroalign

from config import Configuration
from libraries.utils import Utils
from libraries.photometry import Photometry
from astropy.io import fits
import astroalign as aa
import pandas as pd
import numpy as np
import warnings
from PyAstronomy import pyasl
from astropy import time
from dateutil import parser
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


class Lightcurves:

    @staticmethod
    def psf_phot(star):
        """ This function is the wrapper function for the PSF photometry.
        :parameter star - The star for the data to be reduced

        :return star_list - The data frame with the star information
        """

        # get the file list for the photometry
        files = Utils.get_file_list(Configuration.CLEAN_DIRECTORY, Configuration.FILE_EXTENSION)

        # make a data frame with the image information
        images = pd.DataFrame(columns=['filepath', 'jd', 'exptime'])
        lc = pd.DataFrame(columns=['jd', 'mag', 'er', 'x', 'y'])
        bad = 0

        # loop through each image
        for idx, file in enumerate(files):
            skip = 0
            # read in the image
            img = fits.getdata(Configuration.CLEAN_DIRECTORY + file)
            header = fits.getheader(Configuration.CLEAN_DIRECTORY + file)
            dt = parser.parse(header['DATE-OBS'])
            tm = time.Time(dt)
            jd = tm.jd

            # get the stars from the hand picked centroids
            if idx == 0:
                xy_list = pd.read_csv(Configuration.CENTROID_DIRECTORY + file.split('.')[0] + '_list.txt',
                                      delimiter=' ', names=['x', 'y'])
                img0 = img
                src = list()
                for idy, row in xy_list.iterrows():
                    ppair = (row.x, row.y)
                    src.append(ppair)

                src = np.array(src)
            else:
                try:

                    img_transf, (i0_list, i1_list) = aa.find_transform(img0, img)
                    img_calc = aa.matrix_transform(src, img_transf.params)
                    xy_list['x'] = img_calc[:, 0]
                    xy_list['y'] = img_calc[:, 1]

                except astroalign.MaxIterError:
                    skip = 1

            if skip == 0:
                # get the updated PSF by forcing the position and get the new centroids
                star_list = Photometry.force_position(img, xy_list)

                # do the PSF photometry for the stars
                phot_hold = Photometry.single_frame_psf_photometry(star_list, img,
                                                                   Configuration.CLEAN_DIRECTORY + '/' +
                                                                   file.split('.')[0])
                if idx == 0:
                    phot_hold['jd'] = jd
                    phot_hold['hjd'] = pyasl.helio_jd(jd - 2.4e6, Configuration.RA, Configuration.DEC, ) + 2.4e6
                    ph = ((jd - Configuration.TC) / Configuration.PERIOD) % 1
                    if ph > 0.5:
                        ph = ph - 1
                    phot_hold['ph'] = ph
                    phot_chk = phot_hold.copy().reset_index(drop=True)

                else:
                    phot_hold['jd'] = jd
                    phot_hold['hjd'] = pyasl.helio_jd(jd - 2.4e6, Configuration.RA, Configuration.DEC, ) + 2.4e6
                    ph = ((jd - Configuration.TC) / Configuration.PERIOD) % 1
                    if ph > 0.5:
                        ph = ph - 1
                    phot_hold['ph'] = ph
                    phot_chk = phot_chk.append(phot_hold).reset_index(drop=True)

            else:
                bad = bad + 1

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
        files = Utils.get_file_list(Configuration.CLEAN_DIRECTORY, Configuration.FILE_EXTENSION)

        bad = 0

        # loop through each image
        for idx, file in enumerate(files):
            skip = 0
            # read in the image
            img = fits.getdata(Configuration.CLEAN_DIRECTORY + file)
            header = fits.getheader(Configuration.CLEAN_DIRECTORY + file)
            dt = parser.parse(header['DATE-OBS'])
            tm = time.Time(dt)
            jd = tm.jd

            # get the stars from the hand picked centroids
            if idx == 0:
                xy_list = pd.read_csv(Configuration.CENTROID_DIRECTORY + file.split('.')[0] + '_list.txt',
                                      delimiter=' ', names=['x', 'y'])
                img0 = img
                src = list()
                for idy, row in xy_list.iterrows():
                    ppair = (row.x, row.y)
                    src.append(ppair)

                src = np.array(src)
            else:
                try:

                    img_transf, (i0_list, i1_list) = aa.find_transform(img0, img)
                    img_calc = aa.matrix_transform(src, img_transf.params)
                    xy_list['x'] = img_calc[:, 0]
                    xy_list['y'] = img_calc[:, 1]

                except astroalign.MaxIterError:
                    skip = 1

            if skip == 0:
                # get the updated PSF by forcing the position and get the new centroids
                star_list = Photometry.force_position(img, xy_list)

                # do the aperture photometry for the stars
                phot_hold = Photometry.single_frame_aperture_photometry(star_list, img,
                                                                        Configuration.CLEAN_DIRECTORY + '/' +
                                                                        file.split('.')[0], bkg_sub='local')
                if idx == 0:
                    phot_hold['jd'] = jd
                    phot_hold['hjd'] = pyasl.helio_jd(jd - 2.4e6, Configuration.RA, Configuration.DEC, ) + 2.4e6
                    ph = ((jd - Configuration.TC) / Configuration.PERIOD) % 1
                    if ph > 0.5:
                        ph = ph - 1
                    phot_hold['ph'] = ph
                    phot_chk = phot_hold.copy().reset_index(drop=True)
                else:
                    phot_hold['jd'] = jd
                    phot_hold['hjd'] = pyasl.helio_jd(jd - 2.4e6, Configuration.RA, Configuration.DEC, ) + 2.4e6
                    ph = ((jd - Configuration.TC) / Configuration.PERIOD) % 1
                    if ph > 0.5:
                        ph = ph - 1
                    phot_hold['ph'] = ph
                    phot_chk = phot_chk.append(phot_hold).reset_index(drop=True)
            else:
                bad = bad + 1

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
        files = Utils.get_file_list(flux_dir, Configuration.FILE_EXTENSION)

        # combine the flux from the flux files, and write the raw light curves
        Photometry.combine_flux_files(flux_dir, files)

        return

    @staticmethod
    def color_comparison(lc_directory):
        """ This function will compare the colors in a given light curve, and detrend, based on Limbach+2022

        :parameter lc_directory - The directory with the given light curves

        :return A light curve data frame is returned with the updated comparison colors
        """

        # read in the light curve
        lc = pd.read_csv(lc_directory + Configuration.STAR + '_00.lc', header=0, sep=' ')
        tr = pd.read_csv(lc_directory + Configuration.STAR + '_01.lc', header=0, sep=' ')
        lc['cmi_a1'] = lc.apply(lambda x: x.apr1 / np.sum(np.array([x.apr1, x.apr2, x.apr3, x.apr4])), axis=1)
        lc['cmi_a2'] = lc.apply(lambda x: x.apr2 / np.sum(np.array([x.apr1, x.apr2, x.apr3, x.apr4])), axis=1)
        lc['cmi_a3'] = lc.apply(lambda x: x.apr3 / np.sum(np.array([x.apr1, x.apr2, x.apr3, x.apr4])), axis=1)
        lc['cmi_a4'] = lc.apply(lambda x: x.apr4 / np.sum(np.array([x.apr1, x.apr2, x.apr3, x.apr4])), axis=1)

        from sklearn.linear_model import LinearRegression
        X = lc[['apr2', 'apr3', 'apr4']]
        Y = lc[['apr1']]

        regr = LinearRegression()
        res = regr.fit(X, Y)
        lc['mod_ap1'] = res.predict(X)
        lc['c_apr1'] = lc['apr1'] / lc['mod_ap1']

        X = lc[['apr1', 'apr3', 'apr4']]
        Y = lc[['apr2']]

        regr = LinearRegression()
        res = regr.fit(X, Y)
        lc['mod_ap2'] = res.predict(X)
        lc['c_apr2'] = lc['apr2'] / lc['mod_ap2']

        X = lc[['apr1', 'apr2', 'apr4']]
        Y = lc[['apr3']]

        regr = LinearRegression()
        res = regr.fit(X, Y)
        lc['mod_ap3'] = res.predict(X)
        lc['c_apr3'] = lc['apr3'] / lc['mod_ap3']

        X = lc[['apr1', 'apr2', 'apr3']]
        Y = lc[['apr4']]

        regr = LinearRegression()
        res = regr.fit(X, Y)
        lc['mod_ap4'] = res.predict(X)
        lc['c_apr4'] = lc['apr4'] / lc['mod_ap4']

        print('hold')
        return lc