""" This set of functions is primarily used for photometery."""

from config import Configuration
from libraries.utils import Utils
from libraries.preprocessing import Preprocessing
from astropy.io import fits
from photutils.centroids import centroid_sources, centroid_2dg, centroid_quadratic, centroid_com, centroid_1dg
from photutils import aperture_photometry
from photutils.detection import find_peaks
from photutils import EllipticalAperture
from photutils.psf import DAOGroup, BasicPSFPhotometry, extract_stars, EPSFBuilder, IntegratedGaussianPRF
from astropy.modeling import models, fitting
from photutils.aperture import ApertureStats
from photutils.background import MMMBackground
from astropy.modeling.fitting import LevMarLSQFitter, TRFLSQFitter, SLSQPLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm, sigma_clipped_stats
from astropy.table import Table
from astropy.nddata import NDData
from astropy import time
from dateutil import parser
import numpy as np
import pandas as pd
import warnings
from photutils.isophote import EllipseGeometry
from photutils.isophote import Ellipse
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
    def generate_ensemble_psf(image_directory, plot_psf='N'):
        """ This function will build an individual PSF for all of the ETSI bandpases for the given star.

        :parameter image_directory - The directory with the images to generate the PSF from
        :parameter plot_psf - if you want the psf to be plotted
        """

        # get the image names from the directory
        files = Utils.get_file_list(image_directory, '.fits')

        # pull all of the images into a list for PSF creation
        # get the master frame
        nddata = list()
        epsfs = list()
        for file in files:
            img = fits.getdata(image_directory + file)
            img = img - np.median(img)
            nddata.append(NDData(data=img))

        # read in the star positions
        star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.BEAM_TYPE + '_star_list.txt',
                                delimiter=' ', names=['x', 'y'])
        # get the number of stars to look for
        num_stars = int(len(star_list) / Configuration.NUM_PSF)
        star_idx = np.empty(num_stars, dtype=int)

        for idx in range(0, Configuration.NUM_PSF):
            # grab out the index values
            for idy in range(num_stars):
                star_idx[idy] = int(idy * Configuration.NUM_PSF + idx)

            # pull out the necessary PSFs
            star_list_clipped = star_list.iloc[star_idx]
            stars_tbl = Table()
            stars_tbl['x'] = star_list_clipped['x'].to_numpy()
            stars_tbl['y'] = star_list_clipped['y'].to_numpy()

            # make a list with all of the stars
            star_tbl = list()
            for idy in range(len(files)):
                star_tbl.append(stars_tbl)

            # extract the stars from the list
            stars = extract_stars(nddata, star_tbl, size=(91, 91))

            # initialize the EPSF builder
            epsf_builder = EPSFBuilder(oversampling=4, maxiters=3, smoothing_kernel='quadratic')
            epsf, fitted_stars = epsf_builder(stars)
            epsfs.append(epsf)

            # plot the psf if need be
            if plot_psf == 'Y':
                import matplotlib.pyplot as plt
                from astropy.visualization import simple_norm

                norm = simple_norm(epsf.data, 'log', percent=99.)
                plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
                plt.colorbar()
                plt.show()

        return epsfs

    @staticmethod
    def generate_psf(image_directory, find_stars='Y', plot_psf='N'):
        """ This function will perform PSF photometry on given images.

        :parameter image_directory- The file directory with the image to use to generate the PSF
        :parameter find_stars - Y/N if you want the program to automatically find the stars for the PSF, or
        use a star list
        :parameter plot_psf - Y/N if you want to see the PSF plotted

        :return the PSF is returned to use for the data reduction
        """

        if Configuration.SKY_SUBTRACT == 'global':
            mst_nme = '_master_no_bkg.fits'
        else:
            mst_nme = '_master.fits'
        # get the master frame
        img = fits.getdata(image_directory + Configuration.STAR + "_" + Configuration.BEAM_TYPE + mst_nme)

        if find_stars == 'Y':
            # get the PSF peaks on the image
            peaks_tbl = find_peaks(img, threshold=100000)

            # mask images too close to the edge of the frame
            hsize = int(Configuration.PSF_CUTOUT / 2)
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
        stars = extract_stars(nddata, stars_tbl, size=(81, 61))

        epsf_builder = EPSFBuilder(maxiters=1, oversampling=1) # , smoothing_kernel='quadratic')
        epsf, fitted_stars = epsf_builder(stars)

        if plot_psf == 'Y':
            import matplotlib.pyplot as plt
            from astropy.visualization import simple_norm

            norm = simple_norm(epsf.data, 'log', percent=99.)
            plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
            plt.colorbar()
            plt.show()
        return epsf

    @staticmethod
    def centroid_positions(img, star_list):
        """ This function will centroid the l

        :parameter img - The image where the finder will look for stars
        :parameter star_list- The list of hand selected positions

        :return star_list - The list of star positions
        """

        # get the centroids based on the initial xy_list
        if int(Configuration.ELIP_APER_B / 2) % 2 == 0:
            bx_b = int(Configuration.ELIP_APER_B / 2) + 1
        else:
            bx_b = int(Configuration.ELIP_APER_B / 2)

        if int(Configuration.ELIP_APER_A / 2) % 2 == 0:
            bx_a = int(Configuration.ELIP_APER_A / 2) + 1
        else:
            bx_a = int(Configuration.ELIP_APER_A / 2)

        try:
            x_cen, y_cen = centroid_sources(img, star_list.x.to_numpy(), star_list.y.to_numpy(),
                                            box_size=[bx_b, bx_a],
                                            centroid_func=centroid_com)
            star_list = pd.DataFrame(data=zip(x_cen, y_cen), columns=['x', 'y'])
            return star_list
        except ValueError:
            return star_list

    @staticmethod
    def psf_photometry(epsfs, img, star_list, num_psf, phot_type='flux'):
        """ This function will return the PSF photometry of a given image using the EPSF from generated from the master
        frame and convert everything to magnitude.

        :parameter epsfs - A list of EPSF objects to use for photometry
        :parameter img - The image to perform the PSF photometry on
        :parameter star_list - The pandas data frame with the location of the star PSFs
        :parameter num_psf - The number of PSFs for the given throughput
        :parameter phot_type - Output is in magnitude or flux, default is flux
        """

        # get the number of stars to look for
        num_stars = int(len(star_list) / Configuration.NUM_PSF)
        star_idx = np.empty(num_stars, dtype=int)

        flx = np.zeros(len(star_list))
        xfit = np.zeros(len(star_list))
        yfit = np.zeros(len(star_list))

        for idx in range(0, Configuration.NUM_PSF):
            for idy in range(num_stars):
                star_idx[idy] = int(idy * Configuration.NUM_PSF + idx)

            # pull out the necessary PSFs
            star_list_clipped = star_list.iloc[star_idx]

            # convert the star list to the necessary format
            pos = Table(names=['x_0', 'y_0'],
                        data=[star_list_clipped['x'].to_numpy(), star_list_clipped['y'].to_numpy()])

            daogroup = DAOGroup(2.0 * Configuration.SIGMA_PSF * gaussian_sigma_to_fwhm)

            # do the PSF photometry
            mmm_bkg = MMMBackground()
            photometry = BasicPSFPhotometry(group_maker=daogroup,
                                            bkg_estimator=mmm_bkg,
                                            psf_model=epsfs[idx],
                                            fitter=LevMarLSQFitter(),
                                            fitshape=(91, 91))
            phot_table = photometry(image=img, init_guesses=pos)

            for idy in range(num_stars):
                flx[int(idy * Configuration.NUM_PSF + idx)] = phot_table['flux_fit'][idy]
                xfit[int(idy * Configuration.NUM_PSF + idx)] = phot_table['x_fit'][idy]
                yfit[int(idy * Configuration.NUM_PSF + idx)] = phot_table['y_fit'][idy]

        # get the expected error on the psf photometry
        flxe = np.sqrt(flx * Configuration.GAIN)

        # now sum the flux for each star and convert to magnitude
        num_stars = len(flx) // num_psf
        sflx = np.zeros(num_stars)
        sflxe = np.zeros(num_stars)

        # get the white light curve
        for ii in range(0, num_stars):
            for jj in range(0, num_psf):
                sflx[ii] = flx[ii] + flx[ii * num_psf + jj]
                sflxe[ii] = sflxe[ii] + flxe[ii * num_psf + jj]

        # get the expected noise in flux units
        sflxe = np.sqrt(sflxe * Configuration.GAIN)

        return flx, flxe, sflx, sflxe, xfit, yfit

    @staticmethod
    def aperture_photometry(img, header, positions, num_psf, phot_type):
        """ This function will perform aperture photometry on a selected target star to show real time photometry
        outputs
        :parameter img - The "cleaned" image that will have basic photometry run on the data
        :parameter header - The image header to be used for the exposure time correction
        :parameter positions - A data frame with the star positions
        :parameter num_psf - The number of PSFs for the given throughput
        :parameter phot_type - The type of photometry you want output (mag or flux)
        :return mag, mag_e, clr - The combined magnitude and errors are returned as well as the differential color

        """

        xfit = positions['x'].to_numpy()
        yfit = positions['y'].to_numpy()

        positions_abv = positions.copy().reset_index(drop=True)
        positions_blw = positions.copy().reset_index(drop=True)
        positions_abv['y'] = positions_abv['y'] + Configuration.SKY_POS_ABV
        positions_blw['y'] = positions_blw['y'] - Configuration.SKY_POS_BLW

        pos_hold = positions.append(positions_abv).reset_index(drop=True)
        position = pos_hold.append(positions_blw).reset_index(drop=True)

        # set up the apertures for the photometry
        aperture = EllipticalAperture(position,
                                      Configuration.ELIP_APER_A,
                                      Configuration.ELIP_APER_B)

        # run the photometry to get the data table
        phot_table = aperture_photometry(img, aperture, method='exact')

        # subtract the sky background to get the stellar flux
        flx_hold = np.array(phot_table['aperture_sum'])
        flx_star = flx_hold[0:int(len(flx_hold) / 3)]
        flx_abv = flx_hold[int(len(flx_hold) / 3):int(len(flx_hold) / 3) * 2]
        flx_blw = flx_hold[int(len(flx_hold) / 3) * 2:]

        # only use the apertures that had reasonable readings
        bkg = np.nanmean([flx_abv, flx_blw], axis=0)
        flx_max = ApertureStats(img, aperture).max
        # flx_mdn = ApertureStats(img, aperture).median

        # flx_abv = flx_hold[int(len(flx_mdn) / 3):int(len(flx_mdn) / 3) * 2]
        # flx_blw = flx_hold[int(len(flx_mdn) / 3) * 2:]
        # bkg2 = np.nanmean([flx_abv, flx_blw], axis=0) * np.pi * Configuration.ELIP_APER_A * Configuration.ELIP_APER_B

        if Configuration.SKY_SUBTRACT == 'local':

            flx = flx_star - bkg

            # calculate the expected photometric error
            flx_e = np.sqrt(np.array(flx_star) * Configuration.GAIN)

        else:
            # subtract the sky background to get the stellar flux
            flx = flx_star

            # calculate the expected photometric error
            flx_e = np.sqrt(np.array(flx_star) * Configuration.GAIN)

        # now sum the flux for each star and convert to magnitude
        num_stars = len(flx) // num_psf
        sum_flx = np.zeros(num_stars)
        sum_flxe = np.zeros(num_stars)

        # get the white light curve
        for ii in range(0, num_stars):
            for jj in range(0, num_psf):
                sum_flx[ii] = sum_flx[ii] + flx[ii * num_psf + jj]
                sum_flxe[ii] = sum_flxe[ii] + flx_e[ii * num_psf + jj]

        # get the expected noise in flux units
        sum_flx_e = np.sqrt(sum_flxe * Configuration.GAIN)

        # no background subtraction of any kind should take place, make sure exposure time is "not" corrected for to
        # avoid issues with the bias level during the co-add
        if Configuration.SKY_SUBTRACT == 'none':
            exp_tme = 1.
        else:
            exp_tme = header['EXP_TIME']

        # convert to magnitude
        if phot_type == 'mag':
            mag = 25. - 2.5 * np.log10(flx / exp_tme)
            mag_e = (np.log(10.) / 2.5) * (flx_e / flx)
            sum_mag = 25. - 2.5 * np.log10(sum_flx / exp_tme)
            sum_mag_e = (np.log(10.) / 2.5) * (sum_flx_e / sum_flx)
        # keep as flux
        else:
            mag = flx / exp_tme
            mag_e = flx_e
            sum_mag = sum_flx / exp_tme
            sum_mag_e = sum_flx_e / sum_flx

        return mag, mag_e, sum_mag, sum_mag_e, xfit, yfit, bkg, flx_max

    @staticmethod
    def single_band_photometry(image_directory, band):
        """ This function will only perform APERTURE photometery on two bands in the given directory based on the
        hand picked apertures.

        :parameter image_directory - the directory with the images required for photometry
        :parameter band - the band 'aper1, aper2, etc)
        """
        Utils.log("Performing PSF and Aperture photometry on " + Configuration.STAR + " images.", "info")
        # get the file list for the photometry
        files = Utils.get_file_list(image_directory, '.fits')
        num_psf = 2

        # loop through each image
        for idx, file in enumerate(files):
            # read in the image and get the observation time
            img_coadd, header_coadd = fits.getdata(image_directory + file, 0, header=True)

            # get basic information from the image
            jd = np.median([header_coadd['COADD_ST'], header_coadd['COADD_ED']])
            airmass = np.median([header_coadd['AIRMASSS'], header_coadd['AIRMASSE']])
            exp_time = header_coadd['EXP_TIME']

            img = img_coadd
            header = header_coadd
            sky_value = np.median(img)
            sky_sigma = np.std(img)

            # get the stars from the handpicked centroids
            if idx == 0:
                # read in the star list from ds9 hand picking
                star_list_start = pd.read_csv(Configuration.MASTER_DIRECTORY +
                                              Configuration.BEAM_TYPE + '_star_list_' + band + '.txt',
                                              delimiter=' ', names=['x', 'y'])

            # get the updated centroids
            if Configuration.CENTROID == 'Y':
                star_list = Photometry.centroid_positions(img, star_list_start)
            else:
                star_list = star_list_start.copy().reset_index(drop=True)

            aper_mag, aper_er, sum_aper_mag, sum_aper_er, aper_x, aper_y, aper_bkg, aper_max = \
                Photometry.aperture_photometry(img, header, star_list, num_psf, phot_type='flux')

            num_stars = len(star_list) / num_psf

            line_header = 'jd exp_time airmass sky sky_sig aper aper_e ' + \
                          'aper_1 aper_e1 x_a1 y_a1 aper_bkg1 aper_max_1 ' + \
                          'aper_2 aper_e2 x_a2 y_a2 aper_bkg2 aper_max_2\n'

            for star in range(0, int(num_stars)):
                if star < 10:
                    star_id = '0' + str(star)
                else:
                    star_id = str(star)
                line_st = str(np.around(jd, decimals=6)) + ' ' + \
                          str(np.around(exp_time, decimals=0)) + ' ' + \
                          str(np.around(airmass, decimals=3)) + ' ' + \
                          str(np.around(sky_value, decimals=3)) + ' ' + \
                          str(np.around(sky_sigma, decimals=3)) + ' ' + \
                          str(np.around(sum_aper_mag[star], decimals=6)) + ' ' + \
                          str(np.around(sum_aper_er[star], decimals=6))

                big_array = [i for tup in zip(np.around(aper_mag[star * num_psf: (star + 1) * num_psf],
                                                        decimals=6),
                                              np.around(aper_er[star * num_psf: (star + 1) * num_psf],
                                                        decimals=6),
                                              np.around(aper_x[star * num_psf: (star + 1) * num_psf],
                                                        decimals=2),
                                              np.around(aper_y[star * num_psf: (star + 1) * num_psf],
                                                        decimals=2),
                                              np.around(aper_bkg[star * num_psf: (star + 1) * num_psf],
                                                        decimals=2),
                                              np.around(aper_max[star * num_psf: (star + 1) * num_psf],
                                                        decimals=2)) for i in tup]
                big_string = [str(j) for j in big_array]
                line_ed = ' '.join(big_string)
                line = line_st + ' ' + line_ed + '\n'

                if idx == 0:
                    Utils.write_txt(Configuration.LIGHTCURVE_BAND_DIRECTORY +
                                    Configuration.STAR + '_' + Configuration.BEAM_TYPE +
                                    '_' + star_id + '_' + band + '.lc',
                                    'w', line_header)
                Utils.write_txt(Configuration.LIGHTCURVE_BAND_DIRECTORY +
                                Configuration.STAR + '_' + Configuration.BEAM_TYPE +
                                '_' + star_id + '_' + band + '.lc',
                                'a', line)
            if (idx % 100 == 0) & (idx > 0):
                Utils.log("Completed photometry for 100 images for star " + Configuration.STAR + ". " +
                          str(len(files) - idx) + " images remain.", "info")
        return

    @staticmethod
    def psf_and_aperture_photometry():
        """ This function is the wrapper function for the PSF and APERTURE photometry on the image
        :parameter epsf - The PSF produced from the master frame

        :return star_list - The data frame with the star information
        """

        Utils.log("Performing PSF and Aperture photometry on " + Configuration.STAR + " images.", "info")
        # get the file list for the photometry
        if Configuration.CLEAN == 'Y':
            img_dir = Configuration.CLEAN_DIRECTORY
        else:
            if Configuration.COADDED == 'Y':
                img_dir = Configuration.COADD_DIRECTORY
            else:
                img_dir = Configuration.RAW_DIRECTORY
        files = Utils.get_file_list(img_dir, '.fits')

        # loop through each image
        for idx, file in enumerate(files):
            # read in the image and get the observation time
            img_coadd, header_coadd = fits.getdata(img_dir + file, 0, header=True)

            # get basic information from the image
            if Configuration.COADDED == 'Y':
                jd = np.mean([header_coadd['COADD_ST'], header_coadd['COADD_ED']])
                airmass = np.mean([header_coadd['AIRMASSS'], header_coadd['AIRMASSE']])
                exp_time = header_coadd['EXP_TIME']
                if Configuration.BEAM_TYPE == 'transmission':
                    humidity = np.mean([header_coadd['HUMIDS'], header_coadd['HUMIDE']])
                    focus = np.mean([header_coadd['FOCUSS'], header_coadd['FOCUSE']])
                    pressure = np.mean([header_coadd['PRESSS'], header_coadd['PRESSE']])
                    temperature = np.mean([header_coadd['TEMPS'], header_coadd['TEMPE']])
                    dome =  0 # np.mean([header_coadd['DOMES'], header_coadd['DOMEE']])
                    zenith = np.mean([header_coadd['ZENITHS'], header_coadd['ZENITHE']])
                    azimuth = np.mean([header_coadd['AZIMS'], header_coadd['AZIME']])
                else:
                    humidity = 0
                    focus = 0
                    pressure = 0
                    temperature = 0
                    dome = 0
                    zenith = 0
                    azimuth = 0
            else:
                jd = time.Time(parser.parse(header_coadd['DATE-OBS'])).jd
                airmass = header_coadd['AIRMASS']
                exp_time = Configuration.EXPOSURE_TIME
                if Configuration.FILE_EXTENSION == '.fits':
                    focus = header_coadd['FOCUS']
                    pressure = header_coadd['HIERARCH PRESSURE_MB']
                    temperature = header_coadd['HIERARCH TEMPERATURE_C']
                    humidity = header_coadd['HIERARCH HUMIDITY_PERCENT']
                    dome = 0 # header_coadd['HIERARCH DOME_HOUR']
                    zenith = header_coadd['HIERARCH ZD_DEGREE']
                    azimuth = header_coadd['HIERARCH AZ_DEGREE']

            if Configuration.SKY_SUBTRACT == 'global':
                img, header = Preprocessing.sky_subtract(img_coadd, header_coadd, sky_write='N')
                sky_value = header['SKY_MEDN']
                sky_sigma = header['SKY_SIG']
            else:
                img = img_coadd
                header = header_coadd
                sky_value = np.median(img)
                sky_sigma = np.std(img)

            # get the stars from the handpicked centroids
            if idx == 0:
                # read in the star list from ds9 hand picking
                og_star_list = pd.read_csv(Configuration.MASTER_DIRECTORY + Configuration.BEAM_TYPE + '_star_list.txt',
                                           delimiter=' ', names=['x', 'y'])
                ref_img, ref_header = fits.getdata(img_dir + files[0], 0, header=True)
                prv_list = og_star_list.copy().reset_index(drop=True)

                if Configuration.PSF_PHOTOMETRY == 'Y':
                    epsfs = Photometry.generate_ensemble_psf(img_dir, plot_psf='N')

            # get the updated PSF by forcing the position and get the new centroids
            if Configuration.ALIGNMENT == 'Y':
                star_list_start = Preprocessing.calibrate_to_start(ref_img, img, og_star_list, prv_list)
            else:
                star_list_start = prv_list.copy().reset_index(drop=True)

            # get the updated centroids
            if Configuration.CENTROID == 'Y':
                star_list = Photometry.centroid_positions(img, star_list_start)
            else:
                star_list = star_list_start.copy().reset_index(drop=True)
            prv_list = star_list.copy().reset_index(drop=True)

            aper_mag, aper_er, sum_aper_mag, sum_aper_er, aper_x, aper_y, aper_bkg, aper_max = \
                Photometry.aperture_photometry(img, header, star_list,
                                               Configuration.NUM_PSF, phot_type='flux')

            if Configuration.PSF_PHOTOMETRY == 'Y':
                psf_mag, psf_er, sum_psf_mag, sum_psf_er, psf_x, psf_y = \
                    Photometry.psf_photometry(epsfs, img, star_list, Configuration.NUM_PSF, phot_type='flux')

            num_stars = len(star_list) / Configuration.NUM_PSF

            if Configuration.BEAM_TYPE == 'transmission':
                if Configuration.PSF_PHOTOMETRY == 'Y':
                    line_header = 'jd exp_time airmass humidity focus pressure temperature dome zenith azimuth ' \
                                  'sky sky_sig aper aper_e psf psf_e ' + \
                                  'aper_1 aper_e1 x_a1 y_a1 aper_bkg1 aper_max1 psf_1 psf_e1 x_p1 y_p1 ' + \
                                  'aper_2 aper_e2 x_a2 y_a2 aper_bkg2 aper_max2 psf_2 psf_e2 x_p2 y_p2 ' + \
                                  'aper_3 aper_e3 x_a3 y_a3 aper_bkg3 aper_max3 psf_3 psf_e3 x_p3 y_p3 ' + \
                                  'aper_4 aper_e4 x_a4 y_a4 aper_bkg4 aper_max4 psf_4 psf_e4 x_p4 y_p4 ' + \
                                  'aper_5 aper_e5 x_a5 y_a5 aper_bkg5 aper_max5 psf_5 psf_e5 x_p5 y_p5 ' + \
                                  'aper_6 aper_e6 x_a6 y_a6 aper_bkg6 aper_max6 psf_6 psf_e6 x_p6 y_p6 ' + \
                                  'aper_7 aper_e7 x_a7 y_a7 aper_bkg7 aper_max7 psf_7 psf_e7 x_p7 y_p7 ' + \
                                  'aper_8 aper_e8 x_a8 y_a8 aper_bkg8 aper_max8 psf_8 psf_e8 x_p8 y_p8\n'
                else:
                    line_header = 'jd exp_time airmass humidity focus pressure temperature dome zenith azimuth ' \
                                  'sky sky_sig aper aper_e ' + \
                                  'aper_1 aper_e1 x_a1 y_a1 aper_bkg1 aper_max1 ' + \
                                  'aper_2 aper_e2 x_a2 y_a2 aper_bkg2 aper_max2 ' + \
                                  'aper_3 aper_e3 x_a3 y_a3 aper_bkg3 aper_max3 ' + \
                                  'aper_4 aper_e4 x_a4 y_a4 aper_bkg4 aper_max4 ' + \
                                  'aper_5 aper_e5 x_a5 y_a5 aper_bkg5 aper_max5 ' + \
                                  'aper_6 aper_e6 x_a6 y_a6 aper_bkg6 aper_max6 ' + \
                                  'aper_7 aper_e7 x_a7 y_a7 aper_bkg7 aper_max7 ' + \
                                  'aper_8 aper_e8 x_a8 y_a8 aper_bkg8 aper_max8\n'
            else:
                if Configuration.PSF_PHOTOMETRY == 'Y':
                    line_header = 'jd exp_time airmass humidity focus pressure temperature dome zenith azimuth ' \
                                  'sky sky_sig aper aper_e psf psf_e ' + \
                                  'aper_1 aper_e1 x_a1 y_a1 aper_bkg1 aper_max1 psf_1 psf_e1 x_p1 y_p1 ' + \
                                  'aper_2 aper_e2 x_a2 y_a2 aper_bkg2 aper_max2 psf_2 psf_e2 x_p2 y_p2 ' + \
                                  'aper_3 aper_e3 x_a3 y_a3 aper_bkg3 aper_max3 psf_3 psf_e3 x_p3 y_p3 ' + \
                                  'aper_4 aper_e4 x_a4 y_a4 aper_bkg4 aper_max4 psf_4 psf_e4 x_p4 y_p4 ' + \
                                  'aper_5 aper_e5 x_a5 y_a5 aper_bkg5 aper_max5 psf_5 psf_e5 x_p5 y_p5 ' + \
                                  'aper_6 aper_e6 x_a6 y_a6 aper_bkg6 aper_max6 psf_6 psf_e6 x_p6 y_p6 ' + \
                                  'aper_7 aper_e7 x_a7 y_a7 aper_bkg7 aper_max7 psf_7 psf_e7 x_p7 y_p7\n'
                else:
                    line_header = 'jd exp_time airmass humidity focus pressure temperature dome zenith azimuth ' \
                                  'sky sky_sig aper aper_e ' + \
                                  'aper_1 aper_e1 x_a1 y_a1 aper_bkg1 aper_max1 ' + \
                                  'aper_2 aper_e2 x_a2 y_a2 aper_bkg2 aper_max2 ' + \
                                  'aper_3 aper_e3 x_a3 y_a3 aper_bkg3 aper_max3 ' + \
                                  'aper_4 aper_e4 x_a4 y_a4 aper_bkg4 aper_max4 ' + \
                                  'aper_5 aper_e5 x_a5 y_a5 aper_bkg5 aper_max5 ' + \
                                  'aper_6 aper_e6 x_a6 y_a6 aper_bkg6 aper_max6 ' + \
                                  'aper_7 aper_e7 x_a7 y_a7 aper_bkg7 aper_max7\n'

            for star in range(0, int(num_stars)):
                if star < 10:
                    star_id = '0' + str(star)

                if Configuration.PSF_PHOTOMETRY == 'Y':
                    line_st = str(np.around(jd, decimals=6)) + ' ' + \
                              str(np.around(exp_time, decimals=2)) + ' ' + \
                              str(np.around(airmass, decimals=3)) + ' ' + \
                              str(np.around(humidity, decimals=3)) + ' ' + \
                              str(np.around(focus, decimals=0)) + ' ' + \
                              str(np.around(pressure, decimals=3)) + ' ' + \
                              str(np.around(temperature, decimals=3)) + ' ' + \
                              str(np.around(dome, decimals=6)) + ' ' + \
                              str(np.around(zenith, decimals=6)) + ' ' + \
                              str(np.around(azimuth, decimals=6)) + ' ' + \
                              str(np.around(sky_value, decimals=3)) + ' ' + \
                              str(np.around(sky_sigma, decimals=3)) + ' ' + \
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
                                                  np.around(aper_bkg[star * Configuration.NUM_PSF:
                                                                     (star + 1) * Configuration.NUM_PSF],
                                                            decimals=2),
                                                  np.around(aper_max[star * Configuration.NUM_PSF:
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
                else:
                    line_st = str(np.around(jd, decimals=6)) + ' ' + \
                              str(np.around(exp_time, decimals=2)) + ' ' + \
                              str(np.around(airmass, decimals=3)) + ' ' + \
                              str(np.around(humidity, decimals=3)) + ' ' + \
                              str(np.around(focus, decimals=0)) + ' ' + \
                              str(np.around(pressure, decimals=3)) + ' ' + \
                              str(np.around(temperature, decimals=3)) + ' ' + \
                              str(np.around(dome, decimals=6)) + ' ' + \
                              str(np.around(zenith, decimals=6)) + ' ' + \
                              str(np.around(azimuth, decimals=6)) + ' ' + \
                              str(np.around(sky_value, decimals=3)) + ' ' + \
                              str(np.around(sky_sigma, decimals=3)) + ' ' + \
                              str(np.around(sum_aper_mag[star], decimals=6)) + ' ' + \
                              str(np.around(sum_aper_er[star], decimals=6))
                    big_array = [kk for tup in zip(np.around(aper_mag[star * Configuration.NUM_PSF:
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
                                                  np.around(aper_bkg[star * Configuration.NUM_PSF:
                                                                     (star + 1) * Configuration.NUM_PSF],
                                                            decimals=2),
                                                  np.around(aper_max[star * Configuration.NUM_PSF:
                                                                     (star + 1) * Configuration.NUM_PSF],
                                                            decimals=2)) for kk in tup]
                big_string = [str(j) for j in big_array]
                line_ed = ' '.join(big_string)
                line = line_st + ' ' + line_ed + '\n'

                if idx == 0:
                    Utils.write_txt(Configuration.LIGHTCURVE_DIRECTORY +
                                    Configuration.STAR + '_' + Configuration.BEAM_TYPE + '_' + star_id +
                                    '_' + Configuration.SKY_SUBTRACT + '.lc',
                                    'w', line_header)
                Utils.write_txt(Configuration.LIGHTCURVE_DIRECTORY +
                                Configuration.STAR + '_' + Configuration.BEAM_TYPE + '_' + star_id +
                                '_' + Configuration.SKY_SUBTRACT + '.lc',
                                'a', line)
            if (idx % 100 == 0) & (idx > 0):
                Utils.log("Completed photometry for 100 images for star " + Configuration.STAR + ". " +
                          str(len(files) - idx) + " images remain.", "info")
        return
