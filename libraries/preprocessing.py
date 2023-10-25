from config import Configuration
from libraries.utils import Utils
import numpy as np
import os
# from FITS_tools.hcongrid import hcongrid
import astropy
import astropy.stats
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
from astropy.io import fits
import scipy
import scipy.ndimage
from scipy.interpolate import Rbf
from scipy.interpolate import griddata
import astroalign as aa


class Preprocessing:

    @staticmethod
    def calibrate_to_start(img0, img, og_list, pv_list):
        """ This function will calibrate the current image in the real time pipeline to the first image

        :parameter img0 - The first image to be observed
        :parameter img - The current image in the set
        :parameter og_list - The list of the OG pixel positions
        :parameter pv_list - The previous list of positions

        :return star_list - The updated star positions based on the first frame's data
        """

        # set up the coordinate list for the new image
        src = list()
        star_list = og_list.copy().reset_index(drop=True)

        for idy, row in star_list.iterrows():
            ppair = (row.x, row.y)
            src.append(ppair)

        src = np.array(src)
        star_list = star_list.copy().reset_index(drop=True)
        # set up the image transformation
        try:
            # get the transformation offset between the frames
            img_transf, (i0_list, i1_list) = aa.find_transform(img0, img, max_control_points=5,
                                                               detection_sigma=50, min_area=200)
            img_calc = aa.matrix_transform(src, img_transf.params)

            # update the star list with the positions
            star_list['x'] = img_calc[:, 0]
            star_list['y'] = img_calc[:, 1]

            return star_list
        except aa.MaxIterError:
            return pv_list
        except ValueError:
            return pv_list

    @staticmethod
    def mk_bias(image_directory, dark="N", combine_type='median'):
        """ This function will make the master bias frame using the provided image list.
        :parameter image_directory - a directory where the images reside for combination
        :parameter dark - if the file being generated is a dark, update the name
        :parameter combine_type - Either coadd, median, or mean depending on how you want to combine the files

        :return - The bias frame is returned and written to the master directory
        """

        if dark == 'Y':
            file_name = 'dark.fits'
        else:
            file_name = 'bias.fits'

        if os.path.isfile(Configuration.MASTER_DIRECTORY + file_name) == 0:

            if combine_type == 'median':
                # get the image list
                image_list = Utils.get_file_list(image_directory, Configuration.FILE_EXTENSION)

                # determine the number of loops we need to move through for each image
                nfiles = len(image_list)

                # here is the 'holder'
                hold_data = np.ndarray(shape=(nfiles, Configuration.AXIS_X, Configuration.AXIS_Y))

                # update the log
                Utils.log("Generating a master bias/dark frame with median combination. There are " + str(nfiles) +
                          " images to combine.",
                          "info")

                # now loop through the images
                for jj in range(0, nfiles):
                    # read in the image directly into the block_hold
                    if dark == 'Y':
                        bias_tmp = fits.getdata(Configuration.DARKS_DIRECTORY + image_list[jj])
                    else:
                        bias_tmp = fits.getdata(Configuration.BIAS_DIRECTORY + image_list[jj])

                    if np.ndim(bias_tmp) > 2:
                        bias_tmp = bias_tmp[0]

                    hold_data[jj] = bias_tmp

                # median the mini-images into one large image
                bias_image = np.median(hold_data, axis=0)
                std_image = np.std(hold_data, axis=0)

                # pull the header information from the first file of the set
                if dark == 'Y':
                    bias_header = fits.getheader(Configuration.DARKS_DIRECTORY + image_list[0])
                else:
                    bias_header = fits.getheader(Configuration.BIAS_DIRECTORY + image_list[0])

                bias_header['BIAS_COMB'] = 'median'
                bias_header['NUM_BIAS'] = nfiles
                bias_header["BIASSUB"] = 'Y'

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + file_name,
                             bias_image, bias_header, overwrite=True)
                fits.writeto(Configuration.MASTER_DIRECTORY + 'std_'+file_name,
                             std_image, bias_header, overwrite=True)
            elif combine_type == 'coadd':

                num_combine = int(Configuration.BIN_TIME / Configuration.EXPOSURE_TIME)

                # get the image list
                image_list = Utils.get_file_list(image_directory, Configuration.FILE_EXTENSION)

                # determine the number of loops we need to move through for each image
                nfiles = len(image_list)

                # update the log
                Utils.log("Generating a master bias/dark frame with median combination. There are " + str(num_combine) +
                          " images to co-add.",
                          "info")

                # now loop through the images
                for jj in range(0, num_combine):
                    # read in the image directly into the block_hold
                    if dark == 'Y':
                        bias_tmp = fits.getdata(Configuration.DARKS_DIRECTORY + image_list[jj])
                    else:
                        bias_tmp = fits.getdata(Configuration.BIAS_DIRECTORY + image_list[jj])

                    if np.ndim(bias_tmp) > 2:
                        bias_tmp = bias_tmp[0]

                    if jj == 0:
                        bias_image = bias_tmp
                    else:
                        bias_image = bias_image + bias_tmp

                # pull the header information from the first file of the set
                if dark == 'Y':
                    bias_header = fits.getheader(Configuration.DARKS_DIRECTORY + image_list[0])
                else:
                    bias_header = fits.getheader(Configuration.BIAS_DIRECTORY + image_list[0])

                bias_header['BIAS_COMB'] = 'coadd'
                bias_header['NUM_BIAS'] = nfiles

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + file_name,
                             bias_image, bias_header, overwrite=True)
            else:
                # get the image list
                image_list = Utils.get_file_list(image_directory, Configuration.FILE_EXTENSION)

                # determine the number of loops we need to move through for each image
                nfiles = len(image_list)

                # update the log
                Utils.log("Generating a master bias/dark frame with mean combination. There are "
                          + str(nfiles) + " images to combine.", "info")

                for kk in range(0, nfiles):
                    # read in the bias frame
                    if dark == 'Y':
                        bias_tmp = fits.getdata(Configuration.DARKS_DIRECTORY + image_list[kk])
                    else:
                        bias_tmp = fits.getdata(Configuration.BIAS_DIRECTORY + image_list[kk])

                    # initialize if it's the first file, otherwise....
                    if kk == 0:
                        bias_image = bias_tmp
                    else:
                        bias_image = bias_image + bias_tmp

                # generate the mean bias file
                bias_image = bias_image / nfiles

                # pull the header information from the first file of the set
                if dark == 'Y':
                    bias_header = fits.getheader(Configuration.DARKS_DIRECTORY + image_list[kk])
                else:
                    bias_header = fits.getheader(Configuration.BIAS_DIRECTORY + image_list[kk])

                bias_header['BIAS_COMB'] = 'mean'
                bias_header['NUM_BIAS'] = nfiles
                bias_header["BIASSUB"] = 'Y'

                # write the image out to the master directory
                fits.writeto(Configuration.MASTER_DIRECTORY + file_name,
                             bias_image, bias_header, overwrite=True)
        else:
            bias_image = fits.getdata(Configuration.MASTER_DIRECTORY + file_name, 0)

        return bias_image

    @staticmethod
    def mk_flat(image_directory):
        """ This function will make the master flat frame using the provided image list.
        :parameter image_directory - a directory where the images reside for combination

        :return - The bias frame is returned and written to the master directory
        """

        if os.path.isfile(Configuration.MASTER_DIRECTORY + "_flat.fits") == 0:

            # read in the bias file
            bias = Preprocessing.mk_bias(Configuration.BIAS_DIRECTORY)

            # get the image list
            image_list = Utils.get_file_list(image_directory, '.fits')

            # here is the 'holder'
            hold_data = np.zeros((Configuration.AXIS_X, Configuration.AXIS_Y))

            Utils.log("Making the flat frame using a mean.", "info")

            # now loop through the images
            for idx, image in enumerate(image_list):
                # read in the flat file
                flat_tmp = fits.getdata(image_directory + "\\" + image, 0)
                if np.ndim(flat_tmp) > 2:
                    flat_tmp = flat_tmp[0]

                # read in the image directly into the block_hold
                hold_data += flat_tmp - bias

            # median the mini-images into one large image
            flat_image = hold_data / float(len(image_list))
            nflat_image = flat_image / np.median(flat_image)

            # pull the header information from the first file of the set
            flat_header = fits.getheader(image_directory + "\\" + image_list[0])
            flat_header['comb_typ'] = 'mean'
            flat_header['norm_pix'] = 'mean'
            flat_header['num_comb'] = len(image_list)
            flat_header['mean_pix'] = np.mean(nflat_image)
            flat_header['std_pix'] = np.std(nflat_image)
            flat_header['max_pix'] = np.max(nflat_image)
            flat_header['min_pix'] = np.min(nflat_image)
            flat_header['mean'] = np.mean(flat_image)
            flat_header['std'] = np.std(flat_image)
            flat_header['max'] = np.max(flat_image)
            flat_header['min'] = np.min(flat_image)

            # write the image out to the master directory
            fits.writeto(Configuration.CALIBRATION_DATE_DIRECTORY + Configuration.DATE + "_flat.fits",
                         nflat_image, flat_header, overwrite=True)
        else:
            nflat_image = fits.getdata(Configuration.CALIBRATION_DATE_DIRECTORY + Configuration.DATE + "_flat.fits", 0)

        return nflat_image

    @staticmethod
    def sky_subtract(img, header, sky_write='N'):
        """  The function breaks the image into bxs x bxs sections to save memeory, then cleans each section
        the sections are then recombined and smoothed to remove the transitions. The residual image is then subtracted
        from the image and the header is updated appropriately.

        :parameter img - The image to be cleaned
        :parameter header - The header object to be updated
        :parameter sky_write - Y/N if you want to write the residual background for de-bugging

        :return img_sub, header - The cleaned image and updated header file
        """

        # use the sampling space to make the appropriate size vectors
        lop = 2 * Configuration.PIX

        # size holder for later
        sze = int((Configuration.AXIS_X / Configuration.PIX) * (Configuration.AXIS_Y / Configuration.PIX) +
                  (Configuration.AXIS_X / Configuration.PIX) + (Configuration.AXIS_Y / Configuration.PIX) + 1)

        # calculate the sky statistics
        # sky_mean, sky_median, sky_sig = astropy.stats.sigma_clipped_stats(img, sigma=2.5)
        sky_median = np.nanmedian(img)
        sky_sig = np.nanstd(img)

        # create holder arrays for good and bad pixels
        x = np.zeros(shape=sze)
        y = np.zeros(shape=sze)
        v = np.zeros(shape=sze)
        s = np.zeros(shape=sze)
        nd = int(0)

        # begin the sampling of the "local" sky value
        for jj in range(0, Configuration.AXIS_X + Configuration.PIX, Configuration.PIX):
            for kk in range(0, Configuration.AXIS_Y + Configuration.PIX, Configuration.PIX):
                il = np.amax([jj - lop, 0])
                ih = np.amin([jj + lop, Configuration.AXIS_X - 1])
                jl = np.amax([kk - lop, 0])
                jh = np.amin([kk + lop, Configuration.AXIS_Y - 1])
                c = img[jl:jh, il:ih]

                # select the median value with clipping
                # lsky_mean, lsky, ssky = astropy.stats.sigma_clipped_stats(c, sigma=2.5)
                lsky_mean = np.nanmean(c)
                lsky = np.nanmedian(c)
                ssky = np.nanstd(c)
                x[nd] = np.amin([jj, Configuration.AXIS_X - 1])  # determine the pixel to input
                y[nd] = np.amin([kk, Configuration.AXIS_Y - 1])  # determine the pixel to input
                v[nd] = lsky  # median sky
                s[nd] = ssky  # sigma sky
                nd = nd + 1

        # now we want to remove any possible values which have bad sky values
        rj = np.argwhere(v <= 0)  # stuff to remove
        kp = np.argwhere(v > 0)  # stuff to keep

        if len(rj) > 0:

            # keep only the good points
            xgood = x[kp]
            ygood = y[kp]
            vgood = v[kp]

            for jj in range(0, len(rj[0])):
                # select the bad point
                xbad = x[rj[jj]]
                ybad = y[rj[jj]]

                # use the distance formula to get the closest points
                rd = np.sqrt((xgood - xbad) ** 2. + (ygood - ybad) ** 2.)

                # sort the radii
                pp = sorted(range(len(rd)), key=lambda k: rd[k])

                # use the closest 10 points to get a median
                vnear = vgood[pp[0:9]]
                ave = np.median(vnear)

                # insert the good value into the array
                v[rj[jj]] = ave

        # now we want to remove any possible values which have bad sigmas
        rj = np.argwhere(s >= 2 * sky_sig)
        kp = np.argwhere(s < 2 * sky_sig)

        if len(rj) > 0:
            # keep only the good points
            xgood = np.array(x[kp])
            ygood = np.array(y[kp])
            vgood = np.array(v[kp])

            for jj in range(0, len(rj)):
                # select the bad point
                xbad = x[rj[jj]]
                ybad = y[rj[jj]]

                # use the distance formula to get the closest points
                rd = np.sqrt((xgood - xbad) ** 2. + (ygood - ybad) ** 2.)

                # sort the radii
                pp = sorted(range(len(rd)), key=lambda k: rd[k])

                # use the closest 10 points to get a median
                vnear = vgood[pp[0:9]]
                ave = np.median(vnear)
                if np.isfinite(ave) == 0:
                    ave = np.median(v[np.isfinite(v)])

                # insert the good value into the array
                v[rj[jj]] = ave

        # set up a meshgrid to interpolate to
        xi = np.linspace(0, Configuration.AXIS_X - 1, Configuration.AXIS_X)
        yi = np.linspace(0, Configuration.AXIS_Y - 1, Configuration.AXIS_Y)
        xx, yy = np.meshgrid(xi, yi)

        # remove any nan of inf values
        if len(v[~np.isfinite(v)]) > 0:
            v[~np.isfinite(v)] = np.median(v[np.isfinite(v)])

        # now we interpolate to the rest of the image with a cubic interpolation
        res = griddata((x, y), v, (xx, yy), method='cubic')

        # subtract the sky gradient and add back the median background
        img_sub = img - res
        img_sub = img_sub + sky_median

        # attempt "flat-field"
        flat = res / np.median(res)
        img_flt = (img / flat) - res + sky_median

        # update the header
        header['sky_medn'] = sky_median
        header['sky_sig'] = sky_sig
        header['sky_sub'] = 'yes'

        # if desired, write out the sky background to the working directory
        if sky_write == 'Y':
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'sky_background.fits', res, overwrite=True)
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img.fits', img, overwrite=True)
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img_sub.fits', img_sub, overwrite=True)
            fits.writeto(Configuration.ANALYSIS_DIRECTORY + 'img_flt.fits', img_flt, overwrite=True)
        return img_sub, header

    @staticmethod
    def bias_subtract(img, header, dark):
        """ This function will subtract a bias frame
        :parameter img - The image to de-bias / de-dark
        :parameter header - The image header file
        :parameter dark - If Y then the dark frame will be subtracted rather than a bias

        :return bias_sub, header - The updated image and header """

        if dark == 'Y':
            file_name = "dark.fits"
        else:
            file_name = 'bias.fits'

        # read in the bias frame and subtract
        bias = fits.getdata(Configuration.MASTER_DIRECTORY + file_name)

        # subtract the bias from the image
        bias_sub = img - bias

        # update the header
        if dark == 'Y':
            header['dark_sub'] = 'Y'
        else:
            header['bias_sub'] = 'Y'

        return bias_sub, header

    @staticmethod
    def flat_divide(img, header):
        """ This function will divide a flat field.
        :parameter img - The image to flatten
        :parameter header - The image header file

        :return flat_div, header - The updated image and header """
        # read in the flat frame
        flat = fits.getdata(Configuration.MASTER_DIRECTORY + 'flat.fits')

        # subtract the bias from the image
        flat_div = img / flat

        # update the header
        header['flat_div'] = 'Y'

        return flat_div, header

    @staticmethod
    def mk_nme(file, difference_image='N', sky_subtract='N', bias_subtract='N',
               flat_divide='N', alignment='N', dark_subtract="N"):
        """ This function will create the appropriate name for the file based on while steps are taken.
        :argument file - The string with the file name
        :argument difference_image - Y/N if image subtraction occured
        :argument sky_subtract - Y/N if sky subtraction was taken
        :argument bias_subtract - Y/N if a bias is subtracted
        :argument flat_divide - Y/N if a flat field is divided
        :argument alignment - Y/N if the image is aligned
        :argument dark_subtract - Y/N if the image was dark subtracted

        :return file_name - A string with the new file name
        """
        # make a new name for the file based on which actions are taken
        nme_hld = file.split('.')

        # if everything is N then the file name is the original filename
        file_name = file

        # update the file name with a 'd' if at the differencing step
        if difference_image == 'Y':
            file_name = nme_hld[0] + 'd' + Configuration.FILE_EXTENSION

        # otherwise...
        if difference_image == 'N':
            # update the name to be appropriate for what was done to the file
            # nothing occurs
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'N') \
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + Configuration.FILE_EXTENSION
            # bias only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_b' + Configuration.FILE_EXTENSION
            # flat
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_f' + Configuration.FILE_EXTENSION
            # align
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_a' + Configuration.FILE_EXTENSION
            # sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_s' + Configuration.FILE_EXTENSION
            # dark subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_k' + Configuration.FILE_EXTENSION

            # bias and flat only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_bf' + Configuration.FILE_EXTENSION
            # bias and align only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_ba' + Configuration.FILE_EXTENSION
            # bias and sky_subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_bs' + Configuration.FILE_EXTENSION
            # bias and flat and align only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_bfa' + Configuration.FILE_EXTENSION

            # flat and align only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_fa' + Configuration.FILE_EXTENSION
            # flat and sky_subtract only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_fs' + Configuration.FILE_EXTENSION
            # flat and align and sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_fsa' + Configuration.FILE_EXTENSION

            # align and sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_sa' + Configuration.FILE_EXTENSION

            # flat and sky subtract and bias only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'N'):
                file_name = nme_hld[0] + '_bfs' + Configuration.FILE_EXTENSION

            # all steps taken
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bkfsa' + Configuration.FILE_EXTENSION

            # dark and bias only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bk' + Configuration.FILE_EXTENSION
            # dark and flat only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'N') \
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_kf' + Configuration.FILE_EXTENSION
            # dark and align only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_ka' + Configuration.FILE_EXTENSION
            # dark and sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_ks' + Configuration.FILE_EXTENSION
            # dark and bias and flat only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bkf' + Configuration.FILE_EXTENSION
            # dark and bias and align only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bka' + Configuration.FILE_EXTENSION
            # dark and bias and sky subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'N') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bks' + Configuration.FILE_EXTENSION
            # dark and bias and flat and align only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'N')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bkfa' + Configuration.FILE_EXTENSION
            # dark and bias and flat and sky subtract only
            if (bias_subtract == 'Y') and (flat_divide == 'Y') and (alignment == 'N') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_bkfs' + Configuration.FILE_EXTENSION
            # dark and flat and align and sky subtract only
            if (bias_subtract == 'N') and (flat_divide == 'Y') and (alignment == 'Y') and (sky_subtract == 'Y')\
                    and (dark_subtract == 'Y'):
                file_name = nme_hld[0] + '_kfas' + Configuration.FILE_EXTENSION

        return file_name
