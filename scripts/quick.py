""" This script will clean the FFI images, it has been significantly reduced to only
show specific steps. The majority of the functional library can be found in the libraries
directory in utils, fits, photometry, image, and preprocessing libraries."""
from libraries.utils import Utils
from config import Configuration
import os
import cv2
import pandas as pd
from astropy.io import fits
import time
import numpy as np
import warnings
from libraries.photometry import Photometry
from photutils import aperture_photometry
from photutils import EllipticalAperture
from photutils import EllipticalAnnulus
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import logging
logging.getLogger('matplotlib.font_manager').disabled = True
from astropy import time
from dateutil import parser

class Quick:
    @staticmethod
    def clicker_game(file_path):
        # Function which displays an image at file_path and returns the coordinates
        # as a list of tuples of all points clicked until any key is pressed

        # Coordinate initialization
        coords = []

        def clicker_event(event, x, y, flags, params):
            # Internal function defining what to do on a click

            # Check for click
            if event == cv2.EVENT_LBUTTONDOWN:
                # Display coords on image
                font = cv2.FONT_HERSHEY_SIMPLEX
                cv2.putText(img_show, str(x) + ',' +
                            str(y), (x, y), font,
                            1, (0, 0, 0), 2)
                cv2.imshow('Click to record points!', img_show)  # Reload image with text

                # Save coords to list
                coords.append((x, y))

        # Reading in image
        img = fits.getdata(file_path)
        header = fits.getheader(file_path)
        dt = parser.parse(header['DATE-OBS'])
        tm = time.Time(dt)
        time_img  = tm.jd

        ## Small preprocessing to make it legible
        # Removing outliers
        std = np.std(img)
        median = np.median(img)
        noise_lim = .5
        img_clean = np.where(img > median + noise_lim * std, median + noise_lim * std, img)
        img_clean = np.where(img_clean < median - noise_lim * std, median - noise_lim * std, img_clean)

        # Moving distribution to positive values if it has negative values
        if np.min(img_clean) < 0:
            img_clean = img_clean + np.abs(np.min(img_clean))

        # converting to range [0,1]
        img_rescale = img_clean / np.max(img_clean)

        # cv2 only really understands 8bit values
        img_show = np.array(img_rescale * 255, dtype=np.uint8)

        ## Clicker game
        cv2.imshow('Click to record points!', img_show)  # Displaying image
        cv2.setMouseCallback('Click to record points!', clicker_event)  # Creating clicker event

        # Exit condition (press any key)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

        return (coords), img, time_img

    @staticmethod
    def real_time_photometry(img, coords):
        """ This function will perform aperture photometry on a selected target star to show real time photometry
        outputs
        :parameter img - The "cleaned" image that will have basic photometry run on the data
        :parameter coords - The coordinates from clicker.py to use for the analysis

        :return mag, mag_e - The combined magnitude and errors are returned

        """

        # convert the positions to a dataframe
        xy_list = pd.DataFrame(data=coords, columns=['x', 'y'])

        positions = Photometry.force_position(img, xy_list)

        for idx in range(0, Configuration.NUM_PSF):

            # pull out the correct position for the first PSF
            star_loc = positions.iloc[:, idx * 2 + 2:idx * 2 + 4].values.tolist()

            # set up the apertures for the photometry
            aperture = EllipticalAperture(star_loc, Configuration.ELIP_APER_A, Configuration.ELIP_APER_B)
            aperture_annulus = EllipticalAnnulus(star_loc,
                                                 a_in=Configuration.ELIP_ANNULI_A0,
                                                 a_out=Configuration.ELIP_ANNULI_AN,
                                                 b_in=Configuration.ELIP_ANNULI_B0,
                                                 b_out=Configuration.ELIP_ANNULI_BN)
            apers = [aperture, aperture_annulus]

            # run the photometry to get the data table
            phot_table = aperture_photometry(img, apers, method='exact')

            # simple sky calculations
            sky = np.median(img)

            # subtract the sky background to get the stellar flux
            if idx == 0:
                flx = np.array(phot_table['aperture_sum_0'] - sky)
            else:
                flx = flx + np.array(phot_table['aperture_sum_0'] - sky)

            # make sure the appropriate PSFs are being kept for the target store for color
            if idx == (Configuration.COLOR_COMP_1 - 1):
                clr1 = np.array(phot_table['aperture_sum_0'] - sky)
            if idx == (Configuration.COLOR_COMP_2 - 1):
                clr2 = np.array(phot_table['aperture_sum_0'] - sky)

        # calculate the expected photometric error
        star_error = np.sqrt(flx * Configuration.GAIN)
        sky_error = np.sqrt(Configuration.NUM_PSF * aperture.area * sky * Configuration.GAIN)

        # combine sky and signal error in quadrature
        flx_e = np.array(np.sqrt(star_error ** 2 + sky_error ** 2))

        # convert to magnitude
        mag = 25. - 2.5 * np.log10(flx)
        mag_e = (np.log(10.) / 2.5) * (flx_e / flx)

        # get the color change in magnitudes
        clr1_mag = 25. - 2.5 * np.log10(clr1)
        clr2_mag = 25. - 2.5 * np.log10(clr2)
        clr = clr1_mag - clr2_mag

        return mag, mag_e, clr
