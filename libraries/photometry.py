""" This set of functions is primarily used for photometery."""
from config import Configuration
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
from scripts.clean import Clean
from astropy.io import fits
from astropy import time
from dateutil import parser
from photutils.centroids import centroid_sources, centroid_2dg
from photutils import aperture_photometry
from photutils.detection import find_peaks, IRAFStarFinder
from photutils import EllipticalAperture, EllipticalAnnulus
from photutils.psf import DAOGroup, BasicPSFPhotometry, extract_stars, EPSFBuilder, IterativelySubtractedPSFPhotometry
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter, LinearLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm, sigma_clipped_stats
from astropy.table import Table
from astropy.nddata import NDData
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)
import matplotlib.pyplot as plt

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
        file_list = Utils.get_file_list(image_directory, '.fits')
        # img = fits.getdata(Configuration.STAR + "_" + Configuration.BEAM_TYPE + '_master.fits')
        img, header, bd_flag = Clean.clean_img(file_list[0], image_directory,
                                               Configuration.SKY_SUBTRACT,
                                               Configuration.BIAS_SUBTRACT,
                                               Configuration.FLAT_DIVIDE,
                                               Configuration.ALIGNMENT,
                                               Configuration.DARK_SUBTRACT)
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
            star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.BEAM_TYPE + '_star_list.txt',
                                    delimiter=' ', names=['x', 'y'])
            stars_tbl = Table()
            stars_tbl['x'] = star_list['x'].to_numpy()
            stars_tbl['y'] = star_list['y'].to_numpy()

        nddata = NDData(data=img)
        stars = extract_stars(nddata, stars_tbl, size=(61, 81))

        epsf_builder = EPSFBuilder(maxiters=1,
                                   progress_bar=True,
                                   oversampling=1,
                                   shape=(61, 81),
                                   recentering_boxsize=(5, 7),
                                   smoothing_kernel='quartic')
        epsf, fitted_stars = epsf_builder(stars)

        #import matplotlib.pyplot as plt
        #from astropy.visualization import simple_norm

        #norm = simple_norm(epsf.data, 'log', percent=99.)
        #plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        #plt.colorbar()
        #plt.show()
        return epsf

    @staticmethod
    def centroid_positions(img, star_list):
        """ This function will centroid the l

        :parameter img - The image where the finder will look for stars
        :parameter star_list- The list of hand selected positions

        :return star_list - The list of star positions
        """

        # get the centroids based on the initial xy_list
        try:
            x_cen, y_cen = centroid_sources(img, star_list.x.to_numpy(), star_list.y.to_numpy(),
                                            box_size=[Configuration.ELIP_APER_B + 2, Configuration.ELIP_APER_A + 4],
                                            centroid_func=centroid_2dg)
            star_list = pd.DataFrame(data=zip(x_cen, y_cen), columns=['x', 'y'])
            return star_list
        except ValueError:
            return star_list


    @staticmethod
    def psf_photometry(epsf, img, header, star_list, num_psf):
        """ This function will return the PSF photometry of a given image using the EPSF from generated from the master
        frame and convert everything to magnitude.

        :parameter epsf - The PSF function from the master frame
        :parameter img - The image to perform the PSF photometry on
        :parameter header - The exposure time header to use for the correction
        :parameter star_list - The pandas data frame with the location of the star PSFs
        :parameter num_psf - The number of PSFs for the given throughput
        """

        # set up the functions to use for the PSF photometry
        daogroup = DAOGroup(2.0 * Configuration.SIGMA_PSF * gaussian_sigma_to_fwhm)
        bkgrms = MADStdBackgroundRMS()
        std = bkgrms(img)
        mmm_bkg = MMMBackground()
        iraffind = IRAFStarFinder(threshold=3.5 * std,
                                  fwhm=Configuration.SIGMA_PSF * gaussian_sigma_to_fwhm,
                                  minsep_fwhm=0.01, roundhi=5.0, roundlo=-5.0,
                                  sharplo=0.0, sharphi=2.0)

        # convert the star list to the necessary format
        pos = Table(names=['x_0', 'y_0'], data=[star_list['x'].to_numpy(), star_list['y'].to_numpy()])

        # do the PSF photometry
        if Configuration.BKG_FULL_REMOVE == 'Y':
            # photometry = IterativelySubtractedPSFPhotometry(group_maker=daogroup,
            #                                                bkg_estimator=None,
            #                                                psf_model=epsf,
            #                                                fitter=LevMarLSQFitter(),
            #                                                finder=iraffind,
            #                                                niters=3,
            #                                                fitshape=(Configuration.PSF_CUTOUT - 1,
            #                                                          Configuration.PSF_CUTOUT - 1))
            photometry = BasicPSFPhotometry(group_maker=daogroup,
                                            bkg_estimator=None,
                                            psf_model=epsf,
                                            fitter=LevMarLSQFitter(),
                                            fitshape=(61, 81))
        else:
            photometry = BasicPSFPhotometry(group_maker=daogroup,
                                            bkg_estimator=mmm_bkg,
                                            psf_model=epsf,
                                            fitter=LevMarLSQFitter(),
                                            fitshape=(Configuration.PSF_CUTOUT - 1, Configuration.PSF_CUTOUT - 1))
        # generate a table of results
        phot_table = photometry(image=img, init_guesses=pos)
        res_img = photometry.get_residual_image()
        #plt.imshow(res_img, cmap='viridis', aspect=1,
        #           interpolation='nearest', origin='lower')
        #plt.title('Residual Image')
        #plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
        #plt.show()
        fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'sky_background.fits', res_img, overwrite=True)
        #fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img.fits', img, overwrite=True)
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
        mag = 25. - 2.5 * np.log10(flx) - 2.5 * np.log10(header['HIERARCH EXPOSURE TIME'])
        mag_e = (np.log(10.) / 2.5) * (flx_e / flx)
        sum_mag = 25. - 2.5 * np.log10(sum_flx) - 2.5 * np.log10(header['HIERARCH EXPOSURE TIME'])
        sum_mag_e = (np.log(10.) / 2.5) * (sum_flx_e / sum_flx)

        return mag, mag_e, sum_mag, sum_mag_e, xfit, yfit

    @staticmethod
    def aperture_photometry(img, header, star_list, num_psf):
        """ This function will perform aperture photometry on a selected target star to show real time photometry
        outputs
        :parameter img - The "cleaned" image that will have basic photometry run on the data
        :parameter header - The image header to be used for the exposure time correction
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
        if Configuration.BKG_FULL_REMOVE == 'Y':
            flx = np.array(phot_table['aperture_sum_0'])
        else:
            flx = np.array(phot_table['aperture_sum_0']) - \
                  ((phot_table['aperture_sum_1'] / aperture_annulus.area) * aperture.area)

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
        mag = 25. - 2.5 * np.log10(flx) - 2.5 * np.log10(header['HIERARCH EXPOSURE TIME'])
        mag_e = (np.log(10.) / 2.5) * (flx_e / flx)
        sum_mag = 25. - 2.5 * np.log10(sum_flx) - 2.5 * np.log10(header['HIERARCH EXPOSURE TIME'])
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
        if Configuration.MASTER_TYPE == 'normal':
            files = Utils.get_file_list(Configuration.RAW_DIRECTORY, Configuration.FILE_EXTENSION)
        else:
            files = Utils.get_file_list(Configuration.COADD_DIRECTORY, Configuration.FILE_EXTENSION)
        # loop through each image
        for idx, file in enumerate(files):
            # read in the image and get the observation time
            if Configuration.MASTER_TYPE == 'normal':
                img_dirs = Configuration.RAW_DIRECTORY
                img, header, bd_flag = Clean.clean_img(file,
                                                       img_dirs,
                                                       Configuration.SKY_SUBTRACT,
                                                       Configuration.BIAS_SUBTRACT,
                                                       Configuration.FLAT_DIVIDE,
                                                       Configuration.ALIGNMENT,
                                                       Configuration.DARK_SUBTRACT)
            else:
                img_dirs = Configuration.COADD_DIRECTORY
                img, header, bd_flag = Clean.clean_img(file, img_dirs,
                                                       'N',
                                                       'N',
                                                       Configuration.FLAT_DIVIDE,
                                                       Configuration.ALIGNMENT,
                                                       Configuration.DARK_SUBTRACT)
            jd = time.Time(parser.parse(header['DATE-OBS'])).jd

            # get the stars from the handpicked centroids
            if idx == 0:
                # read in the star list from ds9 hand picking
                og_star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.BEAM_TYPE + '_star_list.txt',
                                           delimiter=' ', names=['x', 'y'])
                # save the reference image for later
                # ref_img = fits.getdata(Configuration.MASTER_DIRECTORY +
                #                       Configuration.STAR + '_' + Configuration.BEAM_TYPE + '_master.fits')
                if Configuration.MASTER_TYPE == 'normal':
                    ref_img, ref_header, bd_flag = Clean.clean_img(files[0], img_dirs,
                                                                   Configuration.SKY_SUBTRACT,
                                                                   Configuration.BIAS_SUBTRACT,
                                                                   Configuration.FLAT_DIVIDE,
                                                                   Configuration.ALIGNMENT,
                                                                   Configuration.DARK_SUBTRACT)
                else:
                    ref_img, ref_header, bd_flag = Clean.clean_img(files[0], img_dirs,
                                                                   'N',
                                                                   'N',
                                                                   Configuration.FLAT_DIVIDE,
                                                                   Configuration.ALIGNMENT,
                                                                   Configuration.DARK_SUBTRACT)
                prv_list = og_star_list.copy().reset_index(drop=True)

            # get the updated PSF by forcing the position and get the new centroids
            star_list_start = Preprocessing.calibrate_to_start(ref_img, img, og_star_list, prv_list)

            # get the updated centroids
            star_list = Photometry.centroid_positions(img, star_list_start)
            if (np.sqrt((star_list['x'][0] - prv_list['x'][0]) ** 2 +
                        (star_list['y'][0] - prv_list['y'][0]) ** 2) > 100) & \
                    (np.sqrt((star_list['x'][0] - og_star_list['x'][0]) ** 2 +
                             (star_list['y'][0] - og_star_list['y'][0]) ** 2) > 100) & (idx > 0):

                star_list = prv_list.copy().reset_index(drop=True)

            prv_list = star_list.copy().reset_index(drop=True)
            aper_mag, aper_er, sum_aper_mag, sum_aper_er, aper_x, aper_y = \
                Photometry.aperture_photometry(img, header, star_list, Configuration.NUM_PSF)

            psf_mag, psf_er, sum_psf_mag, sum_psf_er, psf_x, psf_y = \
                Photometry.psf_photometry(epsf, img, header, star_list, Configuration.NUM_PSF)

            num_stars = len(star_list) / Configuration.NUM_PSF

            if Configuration.BEAM_TYPE == 'transmission':
                line_header = 'jd mag_a mage_a mag_p mage_p ' + \
                              'mag_a1 mage_a1 x_a1 y_a1 mag_p1 mage_p1 x_p1 y_p1 ' + \
                              'mag_a2 mage_a2 x_a2 y_a2 mag_p2 mage_p2 x_p2 y_p2 ' + \
                              'mag_a3 mage_a3 x_a3 y_a3 mag_p3 mage_p3 x_p3 y_p3 ' + \
                              'mag_a4 mage_a4 x_a4 y_a4 mag_p4 mage_p4 x_p4 y_p4 ' + \
                              'mag_a5 mage_a5 x_a5 y_a5 mag_p5 mage_p5 x_p5 y_p5 ' + \
                              'mag_a6 mage_a6 x_a6 y_a6 mag_p6 mage_p6 x_p6 y_p6 ' + \
                              'mag_a7 mage_a7 x_a7 y_a7 mag_p7 mage_p7 x_p7 y_p7 ' + \
                              'mag_a8 mage_a8 x_a8 y_a8 mag_p8 mage_p8 x_p8 y_p8\n'
            else:
                line_header = 'jd mag_a mage_a mag_p mage_p ' + \
                              'mag_a1 mage_a1 x_a1 y_a1 mag_p1 mage_p1 x_p1 y_p1 ' + \
                              'mag_a2 mage_a2 x_a2 y_a2 mag_p2 mage_p2 x_p2 y_p2 ' + \
                              'mag_a3 mage_a3 x_a3 y_a3 mag_p3 mage_p3 x_p3 y_p3 ' + \
                              'mag_a4 mage_a4 x_a4 y_a4 mag_p4 mage_p4 x_p4 y_p4 ' + \
                              'mag_a5 mage_a5 x_a5 y_a5 mag_p5 mage_p5 x_p5 y_p5 ' + \
                              'mag_a6 mage_a6 x_a6 y_a6 mag_p6 mage_p6 x_p6 y_p6 ' + \
                              'mag_a7 mage_a7 x_a7 y_a7 mag_p7 mage_p7 x_p7 y_p7\n'

            for star in range(0, int(num_stars)):
                if star < 10:
                    star_id = '0' + str(star)
                line_st = str(np.around(jd, decimals=6)) + ' ' + \
                          str(np.around(sum_aper_mag[star], decimals=6)) + ' ' + \
                          str(np.around(sum_aper_er[star], decimals=6)) + ' ' + \
                          str(np.around(sum_psf_mag[star], decimals=6)) + ' ' + \
                          str(np.around(sum_psf_er[star], decimals=6))
                big_array = [i for tup in zip(np.around(aper_mag[star * Configuration.NUM_PSF:
                                                                 (star + 1) * Configuration.NUM_PSF],
                                                        decimals=6),
                                              np.around(aper_er[star * Configuration.NUM_PSF:
                                                                (star + 1) * Configuration.NUM_PSF],
                                                        decimals=6),
                                              np.around(aper_x[star * Configuration.NUM_PSF:
                                                               (star + 1) * Configuration.NUM_PSF],
                                                        decimals=2),
                                              np.around(aper_y[star * Configuration.NUM_PSF:
                                                               (star + 1) * Configuration.NUM_PSF],
                                                        decimals=2),
                                              np.around(psf_mag[star * Configuration.NUM_PSF:
                                                                (star + 1) * Configuration.NUM_PSF],
                                                        decimals=6),
                                              np.around(psf_er[star * Configuration.NUM_PSF:
                                                               (star + 1) * Configuration.NUM_PSF],
                                                        decimals=6),
                                              np.around(psf_x[star * Configuration.NUM_PSF:
                                                              (star + 1) * Configuration.NUM_PSF],
                                                        decimals=2),
                                              np.around(psf_y[star * Configuration.NUM_PSF:
                                                              (star + 1) * Configuration.NUM_PSF],
                                                        decimals=2)) for i in tup]
                big_string = [str(j) for j in big_array]
                line_ed = ' '.join(big_string)
                line = line_st + ' ' + line_ed + '\n'

                if idx == 0:
                    Utils.write_txt(Configuration.LIGHTCURVE_DIRECTORY +
                                    Configuration.STAR + '_' + Configuration.BEAM_TYPE + '_' + star_id + '.lc',
                                    'w', line_header)
                Utils.write_txt(Configuration.LIGHTCURVE_DIRECTORY +
                                Configuration.STAR + '_' + Configuration.BEAM_TYPE + '_' + star_id + '.lc',
                                'a', line)
            if (idx % 100 == 0) & (idx > 0):
                Utils.log("Completed photometry for 100 images for star " + Configuration.STAR + ". " +
                          str(len(files) - idx) + " images remain.", "info")
        return
