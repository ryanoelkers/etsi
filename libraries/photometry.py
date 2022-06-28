""" This set of functions is primarily used for photometery."""
import os
from config import Configuration
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
from astropy.io import fits
from PyAstronomy import pyasl
from astropy import time
from dateutil import parser
from photutils.centroids import centroid_sources, centroid_2dg
from photutils import aperture_photometry
from photutils.detection import find_peaks
from photutils import EllipticalAperture, EllipticalAnnulus
from photutils.psf import IntegratedGaussianPRF, DAOGroup, BasicPSFPhotometry, extract_stars, EPSFBuilder
from photutils.background import MMMBackground
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.visualization import simple_norm
from astropy.stats import gaussian_sigma_to_fwhm, sigma_clipped_stats
from astropy.table import Table
from astropy.nddata import NDData
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)


class Photometry:


    @staticmethod
    def set_psf_locations(img, xy_list, psf_offset):
        """ This function will set the full positions of each PSF based on an initial prediciton of the separation and
        a centroid to the image.

        :parameter img - The image where the finder will look for stars
        :parameter xy_list - The initial list of the first psf positions
        :parameter psf_offset - The offsets due to the PSF

        :return star_list - A data frame of star list positions
        """
        x = np.array([pos[0] + st for st in psf_offset for pos in xy_list])
        y = np.array([pos[1] + (st * 0) for st in psf_offset for pos in xy_list])

        # get the centroids based on the initial xy_list
        x_cen, y_cen = centroid_sources(img, x, y,
                                        box_size=[Configuration.ELIP_APER_B + 2, Configuration.ELIP_APER_A + 4],
                                        centroid_func=centroid_2dg)
        star_list = pd.DataFrame(data=zip(x_cen, y_cen), columns=['x', 'y'])

        return star_list

    @staticmethod
    def centroid_positions(img, star_list):
        """ This function will centroid the l

        :parameter img - The image where the finder will look for stars
        :parameter star_list- The list of hand selected positions

        :return star_list - The list of star positions
        """

        # get the centroids based on the initial xy_list
        x_cen, y_cen = centroid_sources(img, star_list.x.to_numpy(), star_list.y.to_numpy(),
                                        box_size=[Configuration.ELIP_APER_B + 2, Configuration.ELIP_APER_A + 4],
                                        centroid_func=centroid_2dg)
        star_list = pd.DataFrame(data=zip(x_cen, y_cen), columns=['x', 'y'])

        return star_list

    @staticmethod
    def generate_psf(image_directory, find_stars='Y'):
        """ This function will perform PSF photometry on given images.

        :parameter image_directory- The file directory with the image to use to generate the PSF
        :parameter find_stars - Y/N if you want the program to automatically find the stars for the PSF, or
        use a star list

        :return the PSF is returned to use for the data reduction
        """

        # get the master frame
        img = fits.getdata(image_directory + Configuration.STAR + '_master.fits')

        if find_stars == 'Y':
            # get the PSF peaks on the image
            peaks_tbl = find_peaks(img, threshold=Configuration.PSF_THRESHOLD)

            # mask images too close to the edge of the frame
            hsize = (Configuration.PSF_CUTOUT - 1) / 2
            x = peaks_tbl['x_peak']
            y = peaks_tbl['y_peak']
            mask = ((x > hsize) & (x < (img.shape[1] - 1 - hsize)) &
                    (y > hsize) & (y < (img.shape[0] - 1 - hsize)))

            # select stars for the positions
            stars_tbl = Table()
            stars_tbl['x'] = x[mask]
            stars_tbl['y'] = y[mask]

            # subtract the sky background
            mean_img, median_img, std_img = sigma_clipped_stats(img, sigma=2.)
            img = img - median_img
        else:
            star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + 'star_list.txt',
                                    delimiter=' ', names=['x', 'y'])
            stars_tbl = Table()
            stars_tbl['x'] = star_list['x'].to_numpy()
            stars_tbl['y'] = star_list['y'].to_numpy()

        nddata = NDData(data=img)
        stars = extract_stars(nddata, stars_tbl, size=Configuration.PSF_CUTOUT)

        epsf_builder = EPSFBuilder(maxiters=3, progress_bar=False, oversampling=1)
        epsf, fitted_stars = epsf_builder(stars)

        # import matplotlib.pyplot as plt
        # from astropy.visualization import simple_norm
        # norm = simple_norm(epsf.data, 'log', percent=99.)
        # plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        # plt.colorbar()
        # plt.show()
        return epsf

    @staticmethod
    def single_frame_aperture_photometry(star_list, img, file_name, flux_only='N', bkg_sub='global'):
        """" This function will perform aperture photometry on a single frame.

        :parameter star_list - The list of stars to perform the photometry
        :parameter img - The image where photometry will be performed
        :parameter file_name - The file name for the output flux file
        :parameter flux_only - Y/N if you only want the flux to be generated by the code
        :parameter bkg_sub - global/local/diff if you want to subtract the local background,
                the median image background, or if the frame is from a differenced image

        :return star_list - The star_list data frame is returned with the flux/mag values and errors included
        """

        # set up column names for later
        exs = ['x1', 'x2', 'x3', 'x4']
        wys = ['y1', 'y2', 'y3', 'y4']
        flx = ['flux1', 'flux2', 'flux3', 'flux4']
        flxe = ['fluxe1', 'fluxe2', 'fluxe3', 'fluxe4']
        mag = ['mag1', 'mag2', 'mag3', 'mag4']
        mage = ['mage1', 'mage2', 'mage3', 'mage4']
        lsky = ['lsky1', 'lsky2', 'lsky3', 'lsky4']

        for st in range(0, Configuration.NUM_PSF):
            # generate the positions of the stars
            positions = star_list[[exs[st], wys[st]]].values.tolist()

            # set up the appropraite aperture sizes
            aperture = EllipticalAperture(positions, Configuration.ELIP_APER_A, Configuration.ELIP_APER_B)
            aperture_annulus = EllipticalAnnulus(positions,
                                                 a_in=Configuration.ELIP_ANNULI_A0,
                                                 a_out=Configuration.ELIP_ANNULI_AN,
                                                 b_in=Configuration.ELIP_ANNULI_B0,
                                                 b_out=Configuration.ELIP_ANNULI_BN)
            apers = [aperture, aperture_annulus]

            # run the photometry to get the data table
            phot_table = aperture_photometry(img, apers, method='exact')

            # extract the sky background for each annuli based on either a global or local subtraction
            star_list[lsky[st]] = phot_table['aperture_sum_1'] / aperture_annulus.area

            if bkg_sub == 'local':
                sky = phot_table['aperture_sum_1'] / aperture_annulus.area
            elif bkg_sub == 'global':
                sky = np.median(img)
            else:
                sky = 0

            # subtract the sky background to get the stellar flux and square root of total flux to get the photometric error
            star_list[flx[st]] = np.array(phot_table['aperture_sum_0'] - sky)

            # calculate the expected photometric error
            star_error = np.sqrt((phot_table['aperture_sum_0'] - sky) * Configuration.GAIN)
            sky_error = np.sqrt(aperture.area * sky * Configuration.GAIN)

            # combine sky and signal error in quadrature
            star_list[flxe[st]] = np.array(np.sqrt(star_error ** 2 + sky_error ** 2))

            if flux_only == 'N':
                # convert to magnitude
                star_list[mag[st]] = 25. - 2.5 * np.log10(star_list[flx[st]].to_numpy())
                star_list[mage[st]] = (np.log(10.) / 2.5) * (star_list[flxe[st]].to_numpy() /
                                                             star_list[flx[st]].to_numpy())

        # output the flux files
        star_list.to_csv(file_name + '_'+ Configuration.PHOTOMETRY + '.flux')

        return star_list[star_list.index == 0]

    @staticmethod
    def combine_flux_files(directory, files):
        """ This function combines all of the flux files in a given directory into a single data frame.

        :parameter directory - The directory where the files are located
        :parameter files - A list of the image files

        :returns data_df - A large data frame with all of the stellar flux information
        """

        star_list = pd.read_csv(Configuration.CENTROID_DIRECTORY + files[0].split('.')[0] + '_list.txt',
                                delimiter=' ', names=['x', 'y'])
        num_rrows = len(star_list)

        # make the holders for the light curves
        jdte = np.zeros(len(files))
        hjdte = np.zeros(len(files))
        phse = np.zeros(len(files))

        psf_flx1 = np.zeros((len(files), num_rrows)) - 99.00
        psf_flx2 = np.zeros((len(files), num_rrows)) - 99.00
        psf_flx3 = np.zeros((len(files), num_rrows)) - 99.00
        psf_flx4 = np.zeros((len(files), num_rrows)) - 99.00

        apr_flx1 = np.zeros((len(files), num_rrows)) - 99.00
        apr_flx2 = np.zeros((len(files), num_rrows)) - 99.00
        apr_flx3 = np.zeros((len(files), num_rrows)) - 99.00
        apr_flx4 = np.zeros((len(files), num_rrows)) - 99.00

        x1 = np.zeros((len(files), num_rrows)) - 99.00
        x2 = np.zeros((len(files), num_rrows)) - 99.00
        x3 = np.zeros((len(files), num_rrows)) - 99.00
        x4 = np.zeros((len(files), num_rrows)) - 99.00

        y1 = np.zeros((len(files), num_rrows)) - 99.00
        y2 = np.zeros((len(files), num_rrows)) - 99.00
        y3 = np.zeros((len(files), num_rrows)) - 99.00
        y4 = np.zeros((len(files), num_rrows)) - 99.00

        for idy, file in enumerate(files):

            if os.path.isfile(directory + file.split('.')[0] + '_PSF.flux') & \
                    os.path.isfile(directory + file.split('.')[0] + '_APER.flux'):

                if (idy % 100) == 0:
                    Utils.log("Read in 100 flux files. " + str(len(files)-idy-1) + " files remain.", "info")
                # get the time and phase information
                header = fits.getheader(Configuration.CLEAN_DIRECTORY + file)
                dt = parser.parse(header['DATE-OBS'])
                tm = time.Time(dt)
                jdte[idy] = tm.jd
                hjdte[idy] = pyasl.helio_jd(jdte[idy] - 2.4e6, Configuration.RA, Configuration.DEC) + 2.4e6
                phse[idy] = ((jdte[idy] - Configuration.TC) / Configuration.PERIOD) % 1
                if phse[idy] > 0.5:
                    phse[idy] = phse[idy] - 1

                psf_df = pd.read_csv(directory + file.split('.')[0] + '_PSF.flux', sep=',', index_col=0)
                apr_df = pd.read_csv(directory + file.split('.')[0] + '_APER.flux', sep=',', index_col=0)

                # PSF flux arrays
                psf_flx1[idy, :] = psf_df['flux1'].to_numpy()
                psf_flx2[idy, :] = psf_df['flux2'].to_numpy()
                psf_flx3[idy, :] = psf_df['flux3'].to_numpy()
                psf_flx4[idy, :] = psf_df['flux4'].to_numpy()

                # APER flux arrays
                apr_flx1[idy, :] = apr_df['flux1'].to_numpy()
                apr_flx2[idy, :] = apr_df['flux2'].to_numpy()
                apr_flx3[idy, :] = apr_df['flux3'].to_numpy()
                apr_flx4[idy, :] = apr_df['flux4'].to_numpy()

                # x arrays
                x1[idy, :] = apr_df['x1'].to_numpy()
                x2[idy, :] = apr_df['x2'].to_numpy()
                x3[idy, :] = apr_df['x3'].to_numpy()
                x4[idy, :] = apr_df['x4'].to_numpy()

                # y arrays
                y1[idy, :] = apr_df['y1'].to_numpy()
                y2[idy, :] = apr_df['y2'].to_numpy()
                y3[idy, :] = apr_df['y3'].to_numpy()
                y4[idy, :] = apr_df['y4'].to_numpy()

        # now write hte light curves
        Photometry.write_light_curves(num_rrows, jdte, hjdte, phse,
                                      x1, y1, psf_flx1, apr_flx1,
                                      x2, y2, psf_flx2, apr_flx2,
                                      x3, y3, psf_flx3, apr_flx3,
                                      x4, y4, psf_flx4, apr_flx4)

        return

    @staticmethod
    def write_light_curves(nstars, jd, hjd, ph, x1, y1, pf1, af1, x2, y2, pf2, af2, x3, y3, pf3, af3, x4, y4, pf4, af4):
        """ This function will write the ETSI columns to a text file for later

        :return - Nothing is returned, but the light curve files are written
        """

        # initialize the light curve data frame
        lc = pd.DataFrame(columns={'jd', 'hjd', 'ph',
                                   'x1', 'y1', 'psf1', 'apr1',
                                   'x2', 'y2', 'psf2', 'apr2',
                                   'x3', 'y3', 'psf3', 'apr3',
                                   'x4', 'y4', 'psf4', 'apr4'})

        Utils.log("Starting light curve writing...", "info")

        for idx in range(0, nstars):
            if idx >= 10:
                star_id = str(idx)
            else:
                star_id = '0' + str(idx)

            # add the time, magnitude and error to the data frame
            lc['jd'] = np.around(jd, decimals=6)
            lc['hjd'] = np.around(hjd, decimals=6)
            lc['ph'] = np.around(ph, decimals=6)
            lc['psf1'] = np.around(pf1[:, idx], decimals=6)
            lc['psf2'] = np.around(pf2[:, idx], decimals=6)
            lc['psf3'] = np.around(pf3[:, idx], decimals=6)
            lc['psf4'] = np.around(pf4[:, idx], decimals=6)
            lc['apr1'] = np.around(af1[:, idx], decimals=6)
            lc['apr2'] = np.around(af2[:, idx], decimals=6)
            lc['apr3'] = np.around(af3[:, idx], decimals=6)
            lc['apr4'] = np.around(af4[:, idx], decimals=6)
            lc['x1'] = np.around(x1[:, idx], decimals=6)
            lc['x2'] = np.around(x2[:, idx], decimals=6)
            lc['x3'] = np.around(x3[:, idx], decimals=6)
            lc['x4'] = np.around(x4[:, idx], decimals=6)
            lc['y1'] = np.around(y1[:, idx], decimals=6)
            lc['y2'] = np.around(y2[:, idx], decimals=6)
            lc['y3'] = np.around(y3[:, idx], decimals=6)
            lc['y4'] = np.around(y4[:, idx], decimals=6)

            # write the new file
            lc[['jd', 'hjd', 'ph',
                'x1', 'y1', 'psf1', 'apr1',
                'x2', 'y2', 'psf2', 'apr2',
                'x3', 'y3', 'psf3', 'apr3',
                'x4', 'y4', 'psf4', 'apr4']][lc['apr1'] != -99.000000].to_csv(Configuration.LIGHTCURVE_DIRECTORY +
                                                                              Configuration.STAR + '_' +
                                                                              star_id + ".lc",
                                                                              sep=" ", index=False, na_rep='9.999999')

        return

    @staticmethod
    def centroid_positions(img, star_list):
        """ This function will centroid the l

        :parameter img - The image where the finder will look for stars
        :parameter star_list- The list of hand selected positions

        :return star_list - The list of star positions
        """

        # get the centroids based on the initial xy_list
        x_cen, y_cen = centroid_sources(img, star_list.x.to_numpy(), star_list.y.to_numpy(),
                                        box_size=[Configuration.ELIP_APER_B + 2, Configuration.ELIP_APER_A + 4],
                                        centroid_func=centroid_2dg)
        star_list = pd.DataFrame(data=zip(x_cen, y_cen), columns=['x', 'y'])

        return star_list

    @staticmethod
    def psf_photometry(epsf, img, star_list, num_psf):
        """ This function will return the PSF photometry of a given image using the EPSF from generated from the master
        frame and convert everything to magnitude.

        :parameter epsf - The PSF function from the master frame
        :parameter img - The image to perform the PSF photometry on
        :parameter star_list - The pandas data frame with the location of the star PSFs
        """

        # set up the functions to use for the PSF photometry
        daogroup = DAOGroup(2.0 * Configuration.SIGMA_PSF * gaussian_sigma_to_fwhm)
        mmm_bkg = MMMBackground()

        # convert the star list to the necessary format
        pos = Table(names=['x_0', 'y_0'], data=[star_list['x'].to_numpy(), star_list['y'].to_numpy()])

        # do the PSF photometry
        photometry = BasicPSFPhotometry(group_maker=daogroup,
                                        bkg_estimator=None,
                                        psf_model=epsf,
                                        fitter=LevMarLSQFitter(),
                                        fitshape=(Configuration.PSF_CUTOUT - 1, Configuration.PSF_CUTOUT - 1))

        # generate a table of results
        phot_table = photometry(image=img, init_guesses=pos)

        # pull out the relevant photometry
        flx = np.array(phot_table['flux_fit'])
        try:
            flx_e = np.array(phot_table['flux_unc']) * np.sqrt(Configuration.GAIN)
        except KeyError:
            flx_e = np.array(np.sqrt(phot_table['flux_fit']) * Configuration.GAIN)

        xfit = np.array(phot_table['x_fit'])
        yfit = np.array(phot_table['y_fit'])

        # now sum the flux for each star and convert to magnitude
        num_stars = len(flx) // num_psf
        sum_flx = np.zeros(num_stars)
        sum_flx_e = np.zeros(num_stars)

        # get the white light curve
        for ii in range(0, num_stars):
            for jj in range(0, num_psf):
                sum_flx[ii] = sum_flx[ii] + flx[ii + num_stars * jj]
            sum_flx_e[ii] = np.sqrt(np.sum((flx_e * Configuration.GAIN) ** 2))

        # convert to magnitude
        mag = 25. - 2.5 * np.log10(flx)
        mag_e = (np.log(10.) / 2.5) * (flx_e / flx)
        sum_mag = 25. - 2.5 * np.log10(sum_flx)
        sum_mag_e = (np.log(10.) / 2.5) * (sum_flx_e / sum_flx)

        return mag, mag_e, sum_mag, sum_mag_e, xfit, yfit

    @staticmethod
    def aperture_photometry(img, star_list, num_psf):
        """ This function will perform aperture photometry on a selected target star to show real time photometry
        outputs
        :parameter img - The "cleaned" image that will have basic photometry run on the data
        :parameter star_list - A data frame with the star positions
        :parameter num_psf - The number of PSFs for the given throughput

        :return mag, mag_e, clr - The combined magnitude and errors are returned as well as the differential color

        """

        # convert the positions to a dataframe
        positions = Photometry.centroid_positions(img, star_list)
        xfit = positions['x'].to_numpy()
        yfit = positions['y'].to_numpy()

        # set up the apertures for the photometry
        aperture = EllipticalAperture(positions, Configuration.ELIP_APER_A, Configuration.ELIP_APER_B)
        aperture_annulus = EllipticalAnnulus(positions,
                                             a_in=Configuration.ELIP_ANNULI_A0,
                                             a_out=Configuration.ELIP_ANNULI_AN,
                                             b_in=Configuration.ELIP_ANNULI_B0,
                                             b_out=Configuration.ELIP_ANNULI_BN)
        apers = [aperture, aperture_annulus]

        # run the photometry to get the data table
        phot_table = aperture_photometry(img, apers, method='exact')

        # subtract the sky background to get the stellar flux
        flx = np.array(phot_table['aperture_sum_0'])   # -
                       # (phot_table['aperture_sum_1'] / aperture_annulus.area) * aperture.area)

        # calculate the expected photometric error
        flx_e = np.sqrt(np.array(phot_table['aperture_sum_0']) * Configuration.GAIN)

        # now sum the flux for each star and convert to magnitude
        num_stars = len(flx) // num_psf
        sum_flx = np.zeros(num_stars)
        sum_flxe = np.zeros(num_stars)

        # get the white light curve
        for ii in range(0, num_stars):
            for jj in range(0, num_psf):
                sum_flx[ii] = sum_flx[ii] + flx[ii + num_stars * jj]
                sum_flxe[ii] = sum_flxe[ii] + flx_e[ii + num_stars * jj]

        # get the expected noise in flux units
        sum_flx_e = np.sqrt(sum_flxe * Configuration.GAIN)

        # convert to magnitude
        mag = 25. - 2.5 * np.log10(flx)
        mag_e = (np.log(10.) / 2.5) * (flx_e / flx)
        sum_mag = 25. - 2.5 * np.log10(sum_flx)
        sum_mag_e = (np.log(10.) / 2.5) * (sum_flx_e / sum_flx)

        return mag, mag_e, sum_mag, sum_mag_e, xfit, yfit

    @staticmethod
    def psf_and_aperture_photometry(epsf):
        """ This function is the wrapper function for the PSF and APERTURE photometry on the image
        :parameter epsf - The PSF produced from the master frame

        :return star_list - The data frame with the star information
        """

        Utils.log("Performing PSF and Aperture photometry on " + Configuration.STAR + " images.", "info")
        # get the file list for the photometry
        files = Utils.get_file_list(Configuration.CLEAN_DIRECTORY, Configuration.FILE_EXTENSION)

        # loop through each image
        for idx, file in enumerate(files):
            # read in the image and get the observation time
            img = fits.getdata(Configuration.CLEAN_DIRECTORY + file)
            header = fits.getheader(Configuration.CLEAN_DIRECTORY + file)
            jd = time.Time(parser.parse(header['DATE-OBS'])).jd

            # get the stars from the handpicked centroids
            if idx == 0:
                # read in the star list from ds9 hand picking
                og_star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + 'star_list.txt',
                                           delimiter=' ', names=['x', 'y'])
                # save the reference image for later
                ref_img = img

            else:
                # get the updated PSF by forcing the position and get the new centroids
                og_star_list = Preprocessing.calibrate_to_start(ref_img, img, star_list)

            # get the updated centroids
            star_list = Photometry.centroid_positions(img, og_star_list)

            aper_mag, aper_er, sum_aper_mag, sum_aper_er, aper_x, aper_y = \
                Photometry.aperture_photometry(img, star_list, Configuration.NUM_PSF_TRANSMISSION)

            psf_mag, psf_er, sum_psf_mag, sum_psf_er, psf_x, psf_y = \
                Photometry.psf_photometry(epsf, img, star_list, Configuration.NUM_PSF_TRANSMISSION)

            num_stars = len(star_list) / Configuration.NUM_PSF_TRANSMISSION

            line_header = 'jd mag_a mage_a mag_p mage_p ' +\
                          'mag_a1 mage_a1 x_a1 y_a1 mag_p1 mage_p1 x_p1 y_p1 ' +\
                          'mag_a2 mage_a2 x_a2 y_a2 mag_p2 mage_p2 x_p2 y_p2 '+\
                          'mag_a3 mage_a3 x_a3 y_a3 mag_p3 mage_p3 x_p3 y_p3 '+\
                          'mag_a4 mage_a4 x_a4 y_a4 mag_p4 mage_p4 x_p4 y_p4 '+\
                          'mag_a5 mage_a5 x_a5 y_a5 mag_p5 mage_p5 x_p5 y_p5 '+\
                          'mag_a6 mage_a6 x_a6 y_a6 mag_p6 mage_p6 x_p6 y_p6 '+\
                          'mag_a7 mage_a7 x_a7 y_a7 mag_p7 mage_p7 x_p7 y_p7 '+\
                          'mag_a8 mage_a8 x_a8 y_a8 mag_p8 mage_p8 x_p8 y_p8\n'

            for star in range(0, int(num_stars)):
                if star < 10:
                    star_id = '0' + str(star)
                line_st = str(np.around(jd, decimals=6)) + ' ' + \
                          str(np.around(sum_aper_mag[star], decimals=6)) + ' ' + \
                          str(np.around(sum_aper_er[star], decimals=6)) + ' ' + \
                          str(np.around(sum_psf_mag[star], decimals=6)) + ' ' + \
                          str(np.around(sum_psf_er[star], decimals=6))
                big_array = [i for tup in zip(np.around(aper_mag[star * Configuration.NUM_PSF_TRANSMISSION:
                                                                 (star + 1) * Configuration.NUM_PSF_TRANSMISSION],
                                                        decimals=6),
                                              np.around(aper_er[star * Configuration.NUM_PSF_TRANSMISSION:
                                                                (star + 1) * Configuration.NUM_PSF_TRANSMISSION],
                                                        decimals=6),
                                              np.around(aper_x[star * Configuration.NUM_PSF_TRANSMISSION:
                                                               (star + 1) * Configuration.NUM_PSF_TRANSMISSION],
                                                        decimals=2),
                                              np.around(aper_y[star * Configuration.NUM_PSF_TRANSMISSION:
                                                               (star + 1) * Configuration.NUM_PSF_TRANSMISSION],
                                                        decimals=2),
                                              np.around(psf_mag[star * Configuration.NUM_PSF_TRANSMISSION:
                                                                (star + 1) * Configuration.NUM_PSF_TRANSMISSION],
                                                        decimals=6),
                                              np.around(psf_er[star * Configuration.NUM_PSF_TRANSMISSION:
                                                               (star + 1) * Configuration.NUM_PSF_TRANSMISSION],
                                                        decimals=6),
                                              np.around(psf_x[star * Configuration.NUM_PSF_TRANSMISSION:
                                                              (star + 1) * Configuration.NUM_PSF_TRANSMISSION],
                                                        decimals=2),
                                              np.around(psf_y[star * Configuration.NUM_PSF_TRANSMISSION:
                                                              (star + 1) * Configuration.NUM_PSF_TRANSMISSION],
                                                        decimals=2)) for i in tup]
                big_string = [str(j) for j in big_array]
                line_ed = ' '.join(big_string)
                line = line_st + ' ' + line_ed + '\n'

                if idx == 0:
                    Utils.write_txt(Configuration.LIGHTCURVE_DIRECTORY + 'raw\\' +
                                    Configuration.STAR + '_' + star_id + '.lc',
                                    'w', line_header)
                Utils.write_txt(Configuration.LIGHTCURVE_DIRECTORY + 'raw\\' +
                                Configuration.STAR + '_' + star_id + '.lc',
                                'a', line)
            if (idx % 100 == 0) & (idx > 0):
                Utils.log("Completed photometry for 100 images for star " + Configuration.STAR + ". " +
                          str(len(files) - idx) + " images remain.", "info")
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
            skip = 1
            # read in the image
            img = fits.getdata(Configuration.CLEAN_DIRECTORY + file)
            header = fits.getheader(Configuration.CLEAN_DIRECTORY + file)
            dt = parser.parse(header['DATE-OBS'])
            tm = time.Time(dt)
            jd = tm.jd

            # get the stars from the hand picked centroids
            if idx == 0:
                skip = 0
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

                    img_transf, (i0_list, i1_list) = aa.find_transform(img0, img, min_area=50, max_control_points=10)
                    img_calc = aa.matrix_transform(src, img_transf.params)
                    skip = 0
                    xy_list['x'] = img_calc[:, 0]
                    xy_list['y'] = img_calc[:, 1]

                except astroalign.MaxIterError:
                    print('Bad transform!')
                except ValueError:
                    print('Bad transform!')

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