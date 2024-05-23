""" This set of functions is primarily used for photometery."""
from scripts.calibration import Calibration
from config import Configuration
from libraries.utils import Utils
from astropy.io import fits
from photutils.centroids import centroid_sources, centroid_com
from photutils import aperture_photometry
from photutils import EllipticalAperture
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=Warning)


class Photometry:

    @staticmethod
    def centroid_positions(img, star_list):
        """ This function will centroid the hand-picked PSF positions

        :parameter img - The image where the stars are located
        :parameter star_list- The list of hand selected positions

        :return star_list - The list of centroided star positions
        """

        # get the boxes to centroid at half the expected size for b parameter
        if int(Configuration.ELIP_APER_B / 2) % 2 == 0:
            bx_b = int(Configuration.ELIP_APER_B / 2) + 1
        else:
            bx_b = int(Configuration.ELIP_APER_B / 2)

        # get the boxes to centroid at half the expected size for the a parameter
        if int(Configuration.ELIP_APER_A / 2) % 2 == 0:
            bx_a = int(Configuration.ELIP_APER_A / 2) + 1
        else:
            bx_a = int(Configuration.ELIP_APER_A / 2)

        try:
            # centroid with an updated box size
            x_cen, y_cen = centroid_sources(img, star_list.x.to_numpy(), star_list.y.to_numpy(),
                                            box_size=[bx_b, bx_a],
                                            centroid_func=centroid_com)
            star_list = pd.DataFrame(data=zip(x_cen, y_cen), columns=['x', 'y'])

            return star_list
        # if this fails, return the original star list
        except ValueError:
            return star_list

    @staticmethod
    def aperture_photometry(img, positions):
        """ This function will perform aperture photometry on a selected target star to show real time photometry
        outputs
        :parameter img - The "cleaned" image that will have basic photometry run on the data
        :parameter positions - A data frame with the star positions

        :return flx, flx_e, xfit, yfit, bkg - Return the flux, error, x position, y position, and background
        """

        # get the star positions
        xfit = positions['x'].to_numpy()
        yfit = positions['y'].to_numpy()

        # get the positions of the sky above and below the target positions
        positions_abv = positions.copy().reset_index(drop=True)
        positions_blw = positions.copy().reset_index(drop=True)
        positions_abv['y'] = positions_abv['y'] + Configuration.SKY_POS_ABV
        positions_blw['y'] = positions_blw['y'] - Configuration.SKY_POS_BLW

        # generate a new list of positions
        poss_hold = pd.concat([positions, positions_abv]).reset_index(drop=True)
        position = pd.concat([poss_hold, positions_blw]).reset_index(drop=True)

        # set up the apertures for the photometry
        aperture = EllipticalAperture(position, Configuration.ELIP_APER_A, Configuration.ELIP_APER_B)

        # run the photometry to get the data table
        phot_table = aperture_photometry(img, aperture, method='exact')

        # separate the photometry table into star flux, above bkg flux, and below bkg flux
        flx_hold = np.array(phot_table['aperture_sum'])
        flx_star = flx_hold[0:int(len(flx_hold) / 3)]
        flx_abv = flx_hold[int(len(flx_hold) / 3):int(len(flx_hold) / 3) * 2]
        flx_blw = flx_hold[int(len(flx_hold) / 3) * 2:]

        # only use the apertures that had reasonable readings
        bkg = np.nanmean([flx_abv, flx_blw], axis=0)

        # subtract the sky background to get the stellar flux
        flx = flx_star

        # calculate the expected photon noise from star flux alone
        flx_e = np.sqrt(np.array(flx_star) - np.array(bkg) * Configuration.GAIN)

        return flx, flx_e, xfit, yfit, bkg

    @staticmethod
    def run_photometry(img_dir, wve_str):
        """ This function is the wrapper function for the APERTURE photometry on the co-added images

        :parameter - img_dir - The image directory with the co-added images
        :parameter - wve_str - The name of the wavelengths

        :return - nothing is returned and the light curves are output
        """

        # print out information to the log
        Utils.log("Performing Aperture photometry on " + Configuration.STAR + " images.", "info")

        # get the file list for the photometry
        files = Utils.get_file_list(img_dir, '.fits')

        # loop through each image
        for idx, file in enumerate(files):

            # read in the image
            img_coadd, header_coadd = fits.getdata(img_dir + file, 0, header=True)

            # parse the header for outputs
            jd = np.mean([header_coadd['COADD_ST'], header_coadd['COADD_ED']])
            airmass = np.mean([header_coadd['AIRMASSS'], header_coadd['AIRMASSE']])
            exp_time = header_coadd['EXP_TIME']

            # get global background information
            sky_value = np.median(img_coadd)
            sky_sigma = np.std(img_coadd)

            # get the stars from the handpicked centroids
            if idx == 0:
                # read in the star list from ds9 hand-picking
                og_star_list = pd.read_csv(Configuration.MISC_DIRECTORY +
                                           Configuration.BEAM_TYPE + '_star_list.txt',
                                           delimiter=' ', names=['x', 'y'])

                # this is the first image, so just use the OG star list
                prv_list = og_star_list.copy().reset_index(drop=True)

            # some images have jumps in the data due to movement, if this is the case, align the images to the first one
            if Configuration.ALIGNMENT == 'Y':
                # read in the first image
                ref_img, ref_header = fits.getdata(img_dir + files[0], 0, header=True)

                # align the images to the starting image
                star_list_start, bd_flag = Calibration.calibrate_to_start(ref_img, img_coadd, og_star_list, prv_list)
            else:
                # there is no alignment, so just copy the current positions
                star_list_start = prv_list.copy().reset_index(drop=True)

            # get the updated centroids
            if Configuration.CENTROID == 'Y':
                # centroid the positions
                star_list = Photometry.centroid_positions(img_coadd, star_list_start)
            else:
                # do not centroid the positions
                star_list = star_list_start.copy().reset_index(drop=True)

            # update the previous positions to be the current positions
            prv_list = star_list.copy().reset_index(drop=True)

            # run the aperture photometry
            aper_flx, aper_er, aper_x, aper_y, aper_bkg = Photometry.aperture_photometry(img_coadd, star_list)

            # get the total number of stars in the data
            num_stars = len(star_list) / Configuration.NUM_PSF

            if Configuration.BEAM_TYPE == 'transmission':
                # make the header file for transmission data
                line_header = 'jd exp_time airmass sky sky_sig ' + \
                              wve_str[0] + ' ' + wve_str[0] + '_e ' + wve_str[0] + '_x ' + wve_str[0] + '_y ' + wve_str[0] + '_bkg ' + \
                              wve_str[1] + ' ' + wve_str[1] + '_e ' + wve_str[1] + '_x ' + wve_str[1] + '_y ' + wve_str[1] + '_bkg ' + \
                              wve_str[2] + ' ' + wve_str[2] + '_e ' + wve_str[2] + '_x ' + wve_str[2] + '_y ' + wve_str[2] + '_bkg ' + \
                              wve_str[3] + ' ' + wve_str[3] + '_e ' + wve_str[3] + '_x ' + wve_str[3] + '_y ' + wve_str[3] + '_bkg ' + \
                              wve_str[4] + ' ' + wve_str[4] + '_e ' + wve_str[4] + '_x ' + wve_str[4] + '_y ' + wve_str[4] + '_bkg ' + \
                              wve_str[5] + ' ' + wve_str[5] + '_e ' + wve_str[5] + '_x ' + wve_str[5] + '_y ' + wve_str[5] + '_bkg ' + \
                              wve_str[6] + ' ' + wve_str[6] + '_e ' + wve_str[6] + '_x ' + wve_str[6] + '_y ' + wve_str[6] + '_bkg ' + \
                              wve_str[7] + ' ' + wve_str[7] + '_e ' + wve_str[7] + '_x ' + wve_str[7] + '_y ' + wve_str[7] + '_bkg\n'
            else:
                # make the header file for reflection data
                line_header = 'jd exp_time airmass sky sky_sig ' + \
                              wve_str[0] + ' ' + wve_str[0] + '_e ' + wve_str[0] + '_x ' + wve_str[0] + '_y ' + wve_str[0] + '_bkg ' + \
                              wve_str[1] + ' ' + wve_str[1] + '_e ' + wve_str[1] + '_x ' + wve_str[1] + '_y ' + wve_str[1] + '_bkg ' + \
                              wve_str[2] + ' ' + wve_str[2] + '_e ' + wve_str[2] + '_x ' + wve_str[2] + '_y ' + wve_str[2] + '_bkg ' + \
                              wve_str[3] + ' ' + wve_str[3] + '_e ' + wve_str[3] + '_x ' + wve_str[3] + '_y ' + wve_str[3] + '_bkg ' + \
                              wve_str[4] + ' ' + wve_str[4] + '_e ' + wve_str[4] + '_x ' + wve_str[4] + '_y ' + wve_str[4] + '_bkg ' + \
                              wve_str[5] + ' ' + wve_str[5] + '_e ' + wve_str[5] + '_x ' + wve_str[5] + '_y ' + wve_str[5] + '_bkg ' + \
                              wve_str[6] + ' ' + wve_str[6] + '_e ' + wve_str[6] + '_x ' + wve_str[6] + '_y ' + wve_str[6] + '_bkg\n'

            # move through each star and output the data to text
            for star in range(0, int(num_stars)):

                # generate the star_id (we likely won't have more than 9 stars)
                star_id = '0' + str(star)

                # get the first part of the line based on the image information
                line_st = str(np.around(jd, decimals=6)) + ' ' + \
                          str(np.around(exp_time, decimals=2)) + ' ' + \
                          str(np.around(airmass, decimals=3)) + ' ' + \
                          str(np.around(sky_value, decimals=3)) + ' ' + \
                          str(np.around(sky_sigma, decimals=3))

                # get the star flux information all ready
                big_array = [kk for tup in zip(np.around(aper_flx[star * Configuration.NUM_PSF:
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
                                                         decimals=2)) for kk in tup]

                # make the numerical array a string
                big_string = [str(j) for j in big_array]

                # join the string with a space
                line_ed = ' '.join(big_string)

                # combine the two parts of the line
                line = line_st + ' ' + line_ed + '\n'

                # write out the header if this is the first file
                if idx == 0:
                    Utils.write_txt(Configuration.LIGHTCURVE_DIRECTORY +
                                    Configuration.STAR + '_' + Configuration.BEAM_TYPE + '_' + star_id + '_raw.lc',
                                    'w', line_header)

                # write out the current photometry
                Utils.write_txt(Configuration.LIGHTCURVE_DIRECTORY +
                                Configuration.STAR + '_' + Configuration.BEAM_TYPE + '_' + star_id + '_raw.lc',
                                'a', line)

            # update the terminal as required
            if (idx % 100 == 0) & (idx > 0):
                Utils.log("Completed photometry for 100 images for star " + Configuration.STAR + ". " +
                          str(len(files) - idx) + " images remain.", "info")
        return
